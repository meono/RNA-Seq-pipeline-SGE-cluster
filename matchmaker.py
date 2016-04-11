# A short script that finds Fastq files and pairs them up if they are paired-end reads.
# The script creates a table of the pairs for reference and saves a list of input files in a txt file.
# File type (extension) can be easily modified.

from __future__ import division, print_function
import os
import glob
import tabulate
from tabulate import tabulate

def create_path_if_not_exists(path): # self explanatory
    if not os.path.exists(path):
        os.makedirs(path)


def find_pairs(input_path):

    output_path = os.path.join(input_path, "Fastq_matchmaker")
    create_path_if_not_exists(output_path)

    paired_q = raw_input("Do you have paired reads in separate files? [y/n] ")

    if paired_q.lower() == 'n':
        print('Ok, one file per task.')
        data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))
        # finds any files with fastq.gz extension in folder

    elif paired_q.lower() == 'y':
        print ('Ok, looking for read mates. Mates must be labeled *1.fastq.gz and *2.fastq.gz for pairing.')

        # build forward and reverse file lists, match pairs
        forward = []
        reverse = []
        for read_file in glob.glob(os.path.join(input_path, '*2.fastq.gz')):
            reverse.append(read_file)
        for read_file in glob.glob(os.path.join(input_path, '*1.fastq.gz')):
            forward.append(read_file)

        pairs = []
        for forward_file in forward:
            for_path, ext_for = forward_file.rsplit("_", 1)
            for reverse_file in reverse:
                rev_path, ext_rev = reverse_file.rsplit("_", 1)
                if for_path == rev_path:
                    pairs.append((forward_file, reverse_file))

        # double check correct pairing in  easy to read table
        easy_read = []
        for i in range(len(pairs)):
            path1, name1 = pairs[i][0].rsplit("/", 1)
            path2, name2 = pairs[i][1].rsplit("/", 1)
            easy_read.append((i+1, name1, name2))

        # make a pretty table
        table = tabulate(easy_read, headers=["Pair #", "Forward", "Reverse"], tablefmt="grid")
        print (table)
        mates_q = raw_input("Are these pairings correct? [y/n] ")
        if mates_q.lower() != "y":
            print("Run aborted. Double check file names.")
            return
        else:
            print ("Ok saving table to Matchmaker folder")
            filename = os.path.join(output_path, 'table_pairs.txt')
            f = open(filename, 'w')
            f.write(table)
            f.close()

        data_files = []
        for i in range(len(pairs)):
            data_files.append(str.join(' ', pairs[i]))

    else:
        print("Run aborted.")
        return

    task_count = len(data_files)
    assert task_count > 0, "Could not find any fastq files in folder %s" % input_path


    list_index=os.path.join(output_path,'read_files.txt')
    f=open(list_index,'w')
    f.write(str(data_files))
    f.close()
    print ("All done! Files are located in %s" %output_path)

# get path and run!
input_path = raw_input("Enter the path to your file directory:")
find_pairs(input_path)