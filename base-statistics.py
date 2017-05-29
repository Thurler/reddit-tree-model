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
    filepath = None if out is None else out + "component-hist.png"
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
    # Sort here for plot reasons
    unique = sorted(component_dist.keys())
    counts = [component_dist[x] for x in component_dist]
    total = 1.0 * sum(counts)
    ccdf = map(lambda i: sum(counts[i:]) / total, range(len(counts)))
    # Plot results
    DrawPlot(unique, counts, plt.loglog, "Numero de comentarios", "Frequencia",
             out, out + "comment-component-hist.png")
    DrawPlot(unique, ccdf, plt.loglog, "Numero de comentarios", "CCDF", out,
             out + "comment-component-ccdf.png")


def HeightWidthAnalysis(g, filename, out, force_calc):
    # This should compute height and width sizes for each submission and root
    # comment, and plot the results
    if not force_calc and "width_dist" in g.gp and "height_dist" in g.gp:
        print "Found properties in graph, reading form it!"
        pass  # TODO
    else:
        print "Computing width and height sizes..."

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

    # print "----------------------"

    # print "Analyzing height and width..."
    # HeightWidthAnalysis(g, args.graph, out, args.force_calc)

    print "----------------------"

    print "DONE!"
