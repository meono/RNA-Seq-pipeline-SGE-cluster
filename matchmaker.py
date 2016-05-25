# A short script that finds Fastq files and pairs them up if they are paired-end reads.
# The script creates a table of the pairs for reference and saves a list of input files in a txt file.
# File type (extension) can be easily modified.

from __future__ import division, print_function
import os
import glob
import tabulate
from tabulate import tabulate


def create_path_if_not_exists(mypath):  # self explanatory
    if not os.path.exists(mypath):
        os.makedirs(mypath)


def find_pairs(my_input):

    output_path = os.path.join(my_input, "Fastq_matchmaker")
    create_path_if_not_exists(output_path)

    paired_q = raw_input("Do you have paired reads in separate files? [y/n] ")

    if paired_q.lower() == 'n':
        print('Ok, one file per task.')
        data_files = glob.glob(os.path.join(my_input, '*.fastq.gz'))
        # finds any files with fastq.gz extension in folder

    elif paired_q.lower() == 'y':
        print ('Ok, looking for read mates. Mates must be named *_1.fastq.gz and *_2.fastq.gz for pairing.')

        pairs = []
        full_path = []
        easy_read = []

        for forward_file in glob.glob(os.path.join(my_input, '*_1.fastq.gz')):
            forward_path, forward_name = forward_file.rsplit("/", 1)
            sample_id, ext = forward_name.rsplit("_", 1)
            reverse_name = sample_id + '_2.fastq.gz'

            if os.path.isfile(reverse_name) is True:
                pairs.append((forward_name, reverse_name))
                full_path.append((forward_file, os.path.join(forward_path+"/"+reverse_name)))

        for index, item in enumerate(pairs, start=1):
            easy_read.append((index, item))

        table = tabulate(easy_read, headers=["#", "Pair"], tablefmt="grid")
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

        # correct delimiters for the cluster
        data_files = []
        for i in range(len(full_path)):
            data_files.append(str.join(' ', full_path[i]))

    else:
        print("Run aborted.")
        return

    task_count = len(data_files)
    assert task_count > 0, "Could not find any fastq files in folder %s" % my_input


    list_index=os.path.join(output_path,'read_files.txt')
    f=open(list_index,'w')
    f.write(str(data_files))
    f.close()
    print ("All done! Files are located in %s" %output_path)

# get path and run!
input_path = raw_input("Enter the path to your file directory:")
find_pairs(input_path)