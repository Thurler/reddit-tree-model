#!/usr/bin/env python

import time
import datetime
import os
from os import listdir
from os.path import isfile, join
import argparse
import numpy as np
from graph_tool.all import *

def AddGraph (g1, g2):
    vertex_dic = {}
    for v2 in g2.vertices():
        v1 = g1.add_vertex()
        g1.vp.sid[v1] = g2.vp.sid[v2]
        g1.vp.height[v1] = g2.vp.height[v2]
        vertex_dic[v2] = v1
    for e2 in g2.edges():
        g1.add_edge(vertex_dic[e2.source()], vertex_dic[e2.target()])

if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="Run the R(t,p) model simulator" \
                                                 "for the given arguments.")
    # p and p type
    parser.add_argument("--runs", default=1, type=int,  help="The number of times" \
                        "a model will be generated for each configuration.")
    parser.add_argument("--p_type", default="simple", help="Sets which type of probability" \
                        "will be used. Use 'simple' for a constant p.")
    #
    parser.add_argument("--p_min", type=float, default=0.1, help="Lower bound for p range.")
    parser.add_argument("--p_max", type=float, default=0.9, help="Higher bound for p range.")

    #
    parser.add_argument("--p_step", type=float, default=0.1, help="Step used for incrementing" \
                        "the p value. If p_range = [0.1, 0.9], then a step of 0.1 will generate" \
                        "9 distinct p values.")
    # ttl
    parser.add_argument("--TTL_values", type=list, default=[1000], metavar="[X]", help="TTL configurations.")
    # if true, the iteration pointer gos back to the root when a new node is added, endind the iteration
    parser.add_argument("--add_jump",default=True, help="If set true (default), the iteration" \
                        " pointer goes back to the root when a new node is added, ending the iteration.")
    # Draw?
    parser.add_argument("--draw",default=False, help="If set true (default), draws the graph.")
    parser.add_argument("--debug",default=False, help="If set true (default), enable debug prints.")
    # files
    parser.add_argument("--out_dir", metavar="path/to/dir/", required=True,
                        help="Filename to save resulting graph.")
    args = parser.parse_args()

    # for each p, for each ttl, for each run, execute the simulator once using os.system
    for p in np.arange(args.p_min, args.p_max + args.p_step, args.p_step):
        for ttl in args.TTL_values:
            for i in range(args.runs):
                print ("Running ("+str(i)+") for TTL: "+str(ttl)+ " and p: "+str(p))
                command = "python3 rp-simulator.py --p " + str(p)
                command += " --TTL "+str(ttl)
                command += " --out_dir "+str(args.out_dir)
                if not args.add_jump:
                    command += " --add_jump "+str(args.add_jump)
                if args.debug:
                    command += " --debug "+str(args.debug)
                if args.draw:
                    command += " --draw "+str(args.draw)
                print (command)
                os.system(command)

    print ("Concating .gt files.")
    abs_path = os.path.dirname(os.path.abspath(__file__))
    mypath = os.path.join(abs_path, args.out_dir+"/allgraphs/")
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    # aggregate outputs based on the graph configuration
    config_dic = {}
    for filename in onlyfiles:
        if not filename.endswith(".gt"):
            continue
        cfg = "_".join(filename.split("_")[0:4])
        if cfg not in config_dic:
            config_dic[cfg] = []
        config_dic[cfg].append(filename)

    # add all trees to a single .gt file
    for cfg in config_dic.keys():
        print ("Generating output for "+cfg)
        g = Graph()
        # add properties
        sid_prop = g.new_vertex_property("int16_t")
        height_prop = g.new_vertex_property("int16_t")
        g.vp.sid = sid_prop
        g.vp.height = height_prop
        for filename in config_dic[cfg]:
            sg = load_graph(mypath + filename)
            AddGraph(g, sg)
        out_filename = cfg + "_Final_" + datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S')
        print ("Saving final graph at " + out_filename)
        g.save(os.path.join(abs_path, args.out_dir+"/finalgraphs/" + cfg + ".gt"))


    print ("Batch done.")
