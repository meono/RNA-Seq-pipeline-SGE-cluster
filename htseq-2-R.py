# A script that processes HTSeq-count output for easy input into R.
# Fetches the sample name from the cluster output
# Saves the six lines of basic statistics generated by HTSeq at the bottom of the file in a .txt file

import re
import os
import glob
import argparse
import sys


def create_parser():
    """Return the argument parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_path', required=True,
                        help='''The input folder path, containing the data. ''')
    return parser


def tidy_htseq(my_file, path):
    f = open(my_file, "r")
    for line in f:
        if line.startswith('Input for'):
            fluff, file_name = line.rsplit("/", 1)
            sample, extension = file_name.split(".", 1)

        if re.search(r'\w+\t\d*', line):
            match_stats = re.search(r'^__\w+\t\d*', line)

            if match_stats:
                stats_file = open(os.path.join(path + '/' + sample + "_stats.txt"), "a")
                stats_file.write(line)
                stats_file.close()
            else:
                counts = open(os.path.join(path + '/' + sample + "_htseq-counts.txt"), "a")
                counts.write(line)
                counts.close()
    f.close()
    print ("All done with %s" % sample)


def find_moi(input_path):
    count_files = glob.glob(os.path.join(input_path+"/*.sh.*"))
    print ("Found %d count files" % len(count_files))

    for sample_counts in count_files:
        tidy_htseq(sample_counts, input_path)


def main(argv=None):
    """Program wrapper
    :param argv:
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = create_parser()
    args = parser.parse_args(argv)

    input_path = args.input_path
    find_moi(input_path)

    return 0

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    sys.exit(main())