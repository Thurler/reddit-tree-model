from glob import glob
import argparse
from shutil import copyfile

if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="Extract base statistics from"
                                                 " reddit graph")
    parser.add_argument("--out_dir", metavar="path/to/file", required=True,
                        help="Path where images will be saved")
    parser.add_argument("--in_dir", metavar="path/to/dir", required=True,
                        help="Directory where the images will be searched.")
    parser.add_argument("--plot_filename", required=True,
                        help="The name of the file that will be searched.")
    args = parser.parse_args()

    for dir_path in glob(args.in_dir+"/*/"):
        #print ("DIR "+dir_path)
        for filename in glob(dir_path+"/*"):
            #print ("FILE "+filename)
            if args.plot_filename in filename:
                #print ("TARGET " + filename)

                #get the filename and set the new name as the probability * 100
                new_filename = '%04.f'%(float(dir_path.split("\\")[-2].split("_")[0].replace("simple", ""))*100) + ".png"
                # copy
                copyfile (filename, args.out_dir + "/" + new_filename)
                print ("TARGET "+filename+" COPIED TO "+args.out_dir + "/" + new_filename)

print ("ALL DONE")
