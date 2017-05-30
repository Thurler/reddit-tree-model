#!/usr/bin/env python

import os
import time
import argparse
import numpy as np
from pylab import *
from graph_tool.all import *

DEBUG = False

# Returns the probability function to be used by the simulator.
def ProbabilityTypeFunction (type, base_probability):

    # probability functions are defined here
    def SimpleProbability (node):
        return base_probability;

    if type == 'simple':
        return SimpleProbability
    else:
        raise ValueError('The probability type function name was not recognised.' \
                         'Please check your p_type argument value.')

def ConditionIsMet (TTL, N, vertex_count, iteration_count) :
    if (N < 0) or (N > vertex_count):
        if TTL > iteration_count:
            return True
    return False

def SimulateRtp (p_function, TTL, N, add_jump):
    start_time = time.time()

    # create the graph
    g = Graph()

    # add properties
    sid_prop = g.new_vertex_property("int16_t")
    height_prop = g.new_vertex_property("int16_t")
    g.vp.sid = sid_prop
    g.vp.height = height_prop

    # add the first vertex / root comment
    v_count = 1;
    v1 = g.add_vertex();
    g.vp.sid[v1] = v_count;
    g.vp.height[v1] = 1;

    current_node = v1;
    i_count = 0
    # simulator iterations loop
    while (ConditionIsMet(TTL, N, v_count, i_count)):
        h_count = 1 # height
        # now loop trhough the comment tree
        if (DEBUG):
            print ("=============================")
            print ("Debug at simulator iterations.")
            print ("Current iteration: "+str(i_count))
            print ("Vertex count : "+str(v_count))
        while True:
            if (DEBUG):
                print ("VVV")
                print (current_node)
            replies = []
            degree_sum = 0
            # checks if there will be a reply
            if (p_function(current_node) > np.random.random()):
                # reply! adds new reply (new comment vertex)
                if (DEBUG):
                    print ("Adding new vertex!")
                vn = g.add_vertex()
                g.vp.sid[vn] = v_count
                g.vp.height[vn] = h_count + 1
                g.add_edge(vn, current_node)
                v_count += 1
                # break iteration if add_jump
                if add_jump:
                    current_node = v1
                    break
            # goes to next iteration if there are no replies
            if current_node.in_degree() < 1:
                if (DEBUG):
                    print ("Reached a leaf. Going to next iteration!")
                current_node = v1
                break
            # The next vertex is chosen with probability defined by it's degree,
            # creating a PA relation.
            # First lets select the replies
            if (DEBUG):
                print ("Selecting next vertex.")
                print ("Current vertex indegree: "+ str(current_node.in_degree()))
            for rpl in current_node.in_neighbours():
                # add 1 to give some weight to comments with no replies
                in_degree = rpl.in_degree() + 1
                replies.append((rpl, in_degree))
                degree_sum += in_degree
            # Now lets normalize their degrees
            if (DEBUG):
                # print ("Number of replies: "+ str(len(replies)))
                print ("Degree normaization factor: "+ str(degree_sum))
            replies = [(x[0], float(x[1])/float(degree_sum)) for x in replies]
            # Now choose the reply for the current_node
            p_d = np.random.random()
            p_sum = 0
            for rpl in replies:
                if p_sum + rpl[1] >= p_d:
                    current_node = rpl[0]
                    h_count += 1
                    break
                p_sum += rpl[1]
        i_count += 1

    print ("Simulated Graph Stats:")
    print ("* Number of nodes: " + str(v_count))
    print ("* Number of iterations: " + str(i_count))
    print ("* Time taken: %s seconds" % (time.time() - start_time))
    return g


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="Run the R(t,p) model simulator" \
                                                 "for the given arguments.")
    # p and p type
    parser.add_argument("--p", type=float, metavar=0.5, help="the base probability that the" \
                        "current node will receive a reply", required=True)
    parser.add_argument("--p_type", default="simple", help="sets which type of probability" \
                        "will be used. Use 'simple' for a constant p.")
    # ttl
    parser.add_argument("--TTL", type=int, default=1000, metavar="500", help="Maximum" \
                        "number of iterations. Defaults to 1000.")
    # n
    parser.add_argument("--N", type=int, default=-1, metavar=100, help="Maximum" \
                        "number of nodes to bee added. If not set or negative, will be ignored.")
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

    DEBUG = args.debug

    print ("Starting Simulation.")
    g = SimulateRtp(ProbabilityTypeFunction(args.p_type, args.p), args.TTL, args.N, args.add_jump)

    base_filename = str(args.p) + args.p_type +"_"+ str(args.TTL) +"_"+ str(args.N) +"_"+str(args.add_jump)
    base_filename += "_" + datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S%f')
    abs_path = os.path.dirname(os.path.abspath(__file__))

    print ("Saving graph at " + os.path.join(abs_path, args.out_dir+"/graph/"+base_filename+".gt"))
    g.save(file_name=os.path.join(abs_path, args.out_dir+"/allgraphs/"+base_filename+".gt"))

    if (args.draw):
        print ("Drawing graph.")
        pos = sfdp_layout(g)
        graph_draw(g, pos, output_size=(1000, 1000), vertex_color=[1,1,1,0],
               vertex_size=1, edge_pen_width=1.2,
               vcmap=matplotlib.cm.gist_heat_r, output=abs_path +"/"+args.out_dir+"/plot/"+base_filename+".png")

    print ("All Done")
