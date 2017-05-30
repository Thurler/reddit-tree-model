#!/usr/bin/env python

import argparse
import numpy as np
import graph_tool.all as gt
import matplotlib.pyplot as plt


def LabelComponents(g, filename):
    # Labels out components in graph and saves to file
    components, dist = gt.label_components(g, directed=False)
    gp = g.new_graph_property("vector<int16_t>")
    g.gp["component_dist"] = gp
    g.gp["component_dist"] = dist
    g.vp["conn_components"] = components
    # Save to file so we dont have to compute them again
    print "Saving result in file..."
    g.save(filename)


def DrawPlot(x, y, plotf, xlabel, ylabel, out, filename):
    # Draws a plot based on the plot function and labels it accordingly
    plotf(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # Either output to file or to scree if filename specified
    if out is None:
        plt.show()
    else:
        plt.savefig(filename, bbox_inches='tight')
    # Clear canvas for next plot
    plt.clf()


def PlotDistAnalysis(dist, out, filepath, xlabel):
    unique = sorted(dist.keys())
    counts = [dist[x] for x in unique]
    plus_one = unique[0] == 0
    if plus_one:
        unique = map(lambda x: x+1, unique)
        xlabel += " + 1 (por causa do zero)"
    total = 1.0 * sum(counts)
    ccdf = map(lambda i: sum(counts[i:]) / total, range(len(counts)))
    # Plot results
    filename = None if out is None else out + filepath + "-dist.png"
    DrawPlot(unique, counts, plt.loglog, xlabel, "Frequencia",
             out, filename)
    filename = None if out is None else out + filepath + "-ccdf.png"
    DrawPlot(unique, ccdf, plt.loglog, xlabel, "CCDF", out,
             filename)


def PlotWidthHeightAnalysis(width_dist, submit_width_dist, submit_height_dist,
                            comment_height_dist, out):
    # This should plot all distributions for height and width stuff
    PlotDistAnalysis(submit_width_dist, out, "submits-width", "Largura")
    PlotDistAnalysis(submit_height_dist, out, "submits-height", "Altura")
    PlotDistAnalysis(comment_height_dist, out, "comments-height", "Altura")


def ComponentAnalysis(g, filename, out, force_calc):
    # This should compute and plot how many comments a submission has total
    if not force_calc and "component_dist" in g.gp:
        print "Found property in graph, reading from it!"
    else:
        print "Computing graph components..."
        # Label out components to find connected components (submision trees)
        LabelComponents(g, filename)
    dist = g.gp.component_dist
    # Get distribution of connected component sizes
    unique, counts = np.unique(dist, return_counts=True)
    total = 1.0 * np.sum(counts)
    ccdf = map(lambda i: np.sum(counts[i:]) / total, range(len(counts)))
    # Plot results
    filepath = None if out is None else out + "component-dist.png"
    DrawPlot(unique, counts, plt.loglog, "Numero de comentarios", "Frequencia",
             out, filepath)
    filepath = None if out is None else out + "component-ccdf.png"
    DrawPlot(unique, ccdf, plt.loglog, "Numero de comentarios", "CCDF", out,
             filepath)


def CommentComponentAnalysis(g, filename, out, force_calc):
    # This should compute and plot how many comments each root comment has
    if "conn_components" not in g.vp:
        LabelComponents(g, filename)
    if not force_calc and "comment_component_dist" in g.gp:
        print "Found property in graph, reading from it!"
        component_dist = g.gp.comment_component_dist
    else:
        print "Computing comment components..."
        component_dist = {}
        roots = g.get_in_degrees(g.get_vertices()) == 0
        total_roots = np.sum(roots) - 1
        current = 0
        # We iterate over the vertices that are submissions
        for v in g.get_vertices()[roots]:
            # We generate the subgraph for this submission
            comp_index = g.vp.conn_components[v]
            component = g.new_vertex_property("bool")
            component.a[:] = g.vp.conn_components.a == comp_index
            component.a[v] = 0
            subgraph = gt.GraphView(g, component)
            out_neighbors = g.get_out_neighbours(v)
            total_neighbors = len(out_neighbors)
            print "{} / {} : {}".format(current, total_roots, total_neighbors)
            current += 1
            # For each submission, we iterate over their root comments
            for vv in out_neighbors:
                vv = g.vertex(vv)
                # For each root comment, we count how many vertices there are
                # in their out component (replly tree)
                pmap = gt.label_out_component(subgraph, vv)
                size = int(np.sum(pmap.a))
                # We then archive the distribution for later analysis
                if size not in component_dist:
                    component_dist[size] = 0
                component_dist[size] += 1
        gp = g.new_graph_property("python::object")
        g.gp["comment_component_dist"] = gp
        g.gp["comment_component_dist"] = component_dist
        # Save to file so we dont have to compute this again
        print "Saving result in file..."
        g.save(filename)

    # Plot stuff
    PlotDistAnalysis(component_dist, out, "comment-component",
                     "Numero de comentarios")


def HeightWidthAnalysis(g, filename, out, force_calc):
    # This should compute height and width sizes for each submission and root
    # comment, and plot the results
    if "conn_components" not in g.vp:
        LabelComponents(g, filename)
    if not force_calc and "width_dist" in g.gp:
        print "Found properties in graph, reading form it!"
        width_dist = g.gp.width_dist
        submit_width_dist = g.gp.submit_width_dist
        submit_height_dist = g.gp.submit_height_dist
        comment_height_dist = g.gp.comment_height_dist
    else:
        print "Computing width and height sizes..."
        width_dist = {}
        submit_width_dist = {0: 0}
        submit_height_dist = {0: 0}
        comment_height_dist = {}
        roots = g.get_in_degrees(g.get_vertices()) == 0
        total_roots = np.sum(roots) - 1
        current = 0

        # We use this function to also compute each vertex's level
        vertex_count = len(g.get_vertices())
        vp = g.new_vertex_property("int16_t")

        # We iterate over the vertices that are submissions
        for v in g.get_vertices()[roots]:
            # We generate the subgraph for this submission
            comp_index = g.vp.conn_components[v]
            component = g.new_vertex_property("bool")
            component.a[:] = g.vp.conn_components.a == comp_index
            component.a[v] = 0
            subgraph = gt.GraphView(g, component)
            out_neighbors = g.get_out_neighbours(v)
            total_neighbors = len(out_neighbors)
            print "{} / {} : {}".format(current, total_roots, total_neighbors)
            current += 1
            # Ignore case where there are no neighbors - no width or height
            if total_neighbors == 0:
                submit_width_dist[0] += 1
                submit_height_dist[0] += 1
                continue
            # Store largest height and count width per level
            max_height = 0
            width = {1: total_neighbors}

            # Store this vertex's level as zero
            vp[v] = 0

            # For each submission, we iterate over their root comments
            for vv in out_neighbors:
                # Store this vertex's level as one
                vp[vv] = 1
                vv = g.vertex(vv)
                # For each root comment, we take distance to each other comment
                dist_map = gt.shortest_distance(subgraph, vv, directed=True)
                index_map = np.arange(vertex_count)
                # Filter to those in subgraph
                mask = component.a > 0
                distances = dist_map.a[mask]
                index_map = index_map[mask]
                # Filter to finite distances
                mask = distances < 2147483647
                distances = distances[mask]
                index_map = index_map[mask]
                # Store height information
                height = int(np.max(distances))
                if height not in comment_height_dist:
                    comment_height_dist[height] = 0
                comment_height_dist[height] += 1
                max_height = max(max_height, height + 1)
                # Filter to distances higher than 0
                distances = distances[distances > 0]
                if len(distances) == 0:
                    continue
                # Add one since we excluded root
                distances += 1
                # Store width information
                unique, counts = np.unique(distances, return_counts=True)
                for u, c in zip(unique, counts):
                    if u not in width:
                        width[u] = 0
                    width[u] += c
                # Store level information
                for distance, index in zip(distances, index_map):
                    vp[index] = int(distance)

            # Compute submission height
            if max_height not in submit_height_dist:
                submit_height_dist[max_height] = 0
            submit_height_dist[max_height] += 1
            # Compute submission width
            max_width = sorted(width.items(), key=lambda x: x[1],
                               reverse=True)[0][1]
            if max_width not in submit_width_dist:
                submit_width_dist[max_width] = 0
            submit_width_dist[max_width] += 1
            # Add stuff to width distribution
            for level, w in width.items():
                if level not in width_dist:
                    width_dist[level] = {}
                if w not in width_dist[level]:
                    width_dist[level][w] = 0
                width_dist[level][w] += 1

        # Complete width dist with 0 widths
        for level in width_dist:
            s = 0
            for width in width_dist[level]:
                s += width_dist[level][width]
            width_dist[level][0] = total_roots - s

        # Save level information to file so we don't have to compute it again
        g.vp["levels"] = vp

        # Save distributions to file so we don't have to compute them again
        new_properties = {
            "width_dist": width_dist,
            "submit_width_dist": submit_width_dist,
            "submit_height_dist": submit_height_dist,
            "comment_height_dist": comment_height_dist
        }
        for name, obj in new_properties.items():
            gp = g.new_graph_property("python::object")
            g.gp[name] = gp
            g.gp[name] = obj

        print "Saving result in file..."
        g.save(filename)

    # Plot results
    PlotWidthHeightAnalysis(width_dist, submit_width_dist, submit_height_dist,
                            comment_height_dist, out)


def DegreeAnalysis(g, filename, out, force_calc):
    # This should compute the degree distribution for comments and submissions
    if not force_calc and "submit_deg_dist" in g.gp:
        print "Found properties in graph, reading form it!"
        submit_deg_dist = g.gp.submit_deg_dist
        comment_deg_dist = g.gp.comment_deg_dist
        if "levels" not in g.vp:
            print "Did not find level information for vertices!"
            level_deg_dist = None
        else:
            level_deg_dist = g.gp.level_deg_dist
    else:
        print "Computing degree distributions..."
        submit_deg_dist = {}
        comment_deg_dist = {}
        level_deg_dist = {}
        levels = None if "levels" not in g.vp else g.vp.levels.a
        for vn in g.get_vertices():
            v = g.vertex(vn)
            degree = v.out_degree()
            if v.in_degree() > 0:
                # Not a submission - comment
                if degree not in comment_deg_dist:
                    comment_deg_dist[degree] = 0
                comment_deg_dist[degree] += 1
                # Only do level analysis if level analysis exists
                if levels is not None:
                    level = levels[vn]
                    if level not in level_deg_dist:
                        level_deg_dist[level] = {}
                    if degree not in level_deg_dist[level]:
                        level_deg_dist[level][degree] = 0
                    level_deg_dist[level][degree] += 1
            else:
                # Submission
                if degree not in submit_deg_dist:
                    submit_deg_dist[degree] = 0
                submit_deg_dist[degree] += 1

        # Save distributions to file so we dont have to compute them again
        new_properties = {
            "submit_deg_dist": submit_deg_dist,
            "comment_deg_dist": comment_deg_dist,
            "level_deg_dist": level_deg_dist
        }
        for name, obj in new_properties.items():
            if obj is None:
                continue
            gp = g.new_graph_property("python::object")
            g.gp[name] = gp
            g.gp[name] = obj
        g.save(filename)

    # Plot results
    PlotDistAnalysis(submit_deg_dist, out, "submits-degree",
                     "Numero de comentarios raiz")
    PlotDistAnalysis(comment_deg_dist, out, "comments-degree",
                     "Numero de replies a um comentario")
    for level in level_deg_dist:
        continue
        # something's wrong with levels, which in turn is fucking this up
        # FIX ME
        PlotDistAnalysis(level_deg_dist[level], out,
                         "level-degree-{}".format(level),
                         "Numero de replies no nivel {}".format(level))

if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="Extract base statistics from"
                                                 " reddit graph")
    parser.add_argument("--graph", metavar="path/to/file", required=True,
                        help="Filename to read graph from")
    parser.add_argument("--out", metavar="path/to/dir",
                        help="Directory to save resulting statistics")
    parser.add_argument("--force_calc", default=False, action="store_true",
                        help="Force algorithms to run again")
    args = parser.parse_args()

    print "Loading graph..."
    g = gt.load_graph(args.graph)
    out = None if args.out == "" else args.out
    print "Graph loaded successfully!"

    print "----------------------"

    print "Analyzing submission components..."
    ComponentAnalysis(g, args.graph, out, args.force_calc)

    print "----------------------"

    print "Analyzing comment components..."
    CommentComponentAnalysis(g, args.graph, out, args.force_calc)

    print "----------------------"

    print "Analyzing height and width..."
    HeightWidthAnalysis(g, args.graph, out, args.force_calc)

    print "----------------------"

    print "Analyzing degree distribution..."
    DegreeAnalysis(g, args.graph, out, args.force_calc)

    print "----------------------"

    print "DONE!"
