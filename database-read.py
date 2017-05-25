#!/usr/bin/env python

import argparse
import sqlite3 as sql

from graph_tool.all import *

START_ID = "cqug90g"  # IDs before this one are not in database


def InitGraph(store_comments):
    # Creates vertex properties and assigns them to a new graph
    properties_map = {
        "id": "string",
        "parent_id": "string",
        "score": "int16_t",
        "subreddit": "string"
    }
    if store_comments:
        properties_map["body"] = "string"

    g = Graph()
    for pname, ptype in properties_map.items():
        p = g.new_vertex_property(ptype)
        g.vp[pname] = p

    return g


def ReadDB(filename, subreddit):
    # Read filename as a sqlite database and run reddit query
    conn = sql.connect(filename)
    cursor = conn.cursor()
    if subreddit:
        cursor.execute("SELECT * FROM May2015 WHERE subreddit = ?", subreddit)
    else:
        cursor.execute("SELECT * FROM May2015")
    return cursor


def ParseDB(cursor, graph, store_comments):
    # Parse each db row as a graph vertex
    id_index_map = {}
    submission_id_index_map = {}
    missing_edges = {}
    properties_index_map = {
        "body": 17,
        "id": 9,
        "parent_id": 21,
        "score": 15,
        "subreddit": 8
    }

    i_count = 0
    for row in cursor:
        new_submission = False

        # Each line is a vertex
        graph.add_vertex()

        # Store id reference
        row_id = row[properties_index_map["id"]]
        graph.vp.id[i_count] = row_id
        id_index_map[row_id] = i_count

        # Store score and subreddit info
        subreddit = row[properties_index_map["subreddit"]]
        graph.vp.subreddit[i_count] = subreddit
        graph.vp.score[i_count] = row[properties_index_map["score"]]

        # Store comment body if argument specifies so
        if store_comments:
            graph.vp.body[i_count] = row[properties_index_map["body"]]

        # Check if this id has any missing edges, add them to graph if so
        if row_id in missing_edges:
            for dest in missing_edges[row_id]:
                graph.add_edge(graph.vertex(i_count), graph.vertex(dest))
            del missing_edges[row_id]

        # Store parent id reference
        row_parent_id = row[properties_index_map["parent_id"]]
        prefix, pid = row_parent_id.split('_')

        if prefix == "t3":
            # Parent is a submimssion
            graph.vp.parent_id[i_count] = row_parent_id
            p_index = None
            if pid in submission_id_index_map:
                # Get submission id from map
                p_index = submission_id_index_map[pid]
            else:
                # Submission not in map - add new vertex
                new_submission = True
                graph.add_vertex()
                p_index = i_count + 1
                graph.vp.id[p_index] = row_parent_id
                graph.vp.parent_id[p_index] = ""
                graph.vp.score[p_index] = -1
                graph.vp.subreddit[p_index] = subreddit
                # Add new vertex to submission index map
                submission_id_index_map[pid] = p_index

            # Add edge in graph
            graph.add_edge(graph.vertex(p_index), graph.vertex(i_count))

        elif prefix == "t1" and pid >= START_ID:
            # Parent is a comment within analyzed month
            graph.vp.parent_id[i_count] = pid
            if pid in id_index_map:
                # Get comment id from map and edge to graph
                p_index = id_index_map[pid]
                graph.add_edge(graph.vertex(p_index), graph.vertex(i_count))
            else:
                if pid not in missing_edges:
                    missing_edges[pid] = []
                missing_edges[pid].append(i_count)

        elif pid >= START_ID:
            # Parent is something weird
            print "Error parsing parent id! Row:", row
            return None

        # Increment index counter for next line - extra if new submission
        i_count += (1 + new_submission)

    return len(id_index_map), len(submission_id_index_map)

if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="Read sqlite file as a graph")
    parser.add_argument("--db", metavar="path/to/file", help="Filename to use",
                        required=True)
    parser.add_argument("--out", metavar="path/to/file", required=True,
                        help="Filename to save resulting graph")
    parser.add_argument("--sub", metavar="brasil", help="Subreddit to query")
    parser.add_argument("--store_comments", default=False, action="store_true",
                        help="Store original plaintext comment")
    args = parser.parse_args()

    # Init Graph structure and sqlite cursor
    g = InitGraph(args.store_comments)
    cursor = ReadDB(args.db, args.sub)

    # Parse rows
    comments, submissions = ParseDB(cursor, g, args.store_comments)

    print "LOADED GRAPH SUCCESSFULLY"
    print "{} VERTICES, OF WHICH:".format(comments + submissions)
    print "{} COMMENTS".format(comments)
    print "{} SUBMISSIONS".format(submissions)

    g.save(args.out)
