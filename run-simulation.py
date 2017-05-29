import os
import subprocess
import argparse
import numpy as np

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

    print ("All done.")
