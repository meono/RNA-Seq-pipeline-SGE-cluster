from __future__ import division, print_function
import collections
import logging

logger = logging.getLogger(__name__)

def check_read_uniqueness(reads):
    read_files = [read for read_list in reads.values() for read in read_list]
    if len(read_files) == len(set(read_files)):
        return True
    else:
        report = {}
        duplicates = [read for read, count in collections.Counter(read_files).items() if count > 1]
        for duplicate in duplicates:
            for path, read_list in reads.items():
                if duplicate in read_list:
                    try:
                        report[duplicate].add(path)
                    except:
                        report[duplicate] = {path}
        logger.warning('There are reads with identical names under different paths:\n')
        for read, paths in report.items():
            logger.warning('{} found under'.format(read))
            for path in paths:
                logger.warning('\t{}'.format(path))
        logger.error('Check input reads and restart.')
        return False


def lcs(file1, file2):
    # find longest common substring from *beginning f the string*
    for x in range(1, 1 + len(file1)):
        if file1[0:x] == file2[0:x]:
            tmp = file1[0:x]
    if 'tmp' in locals():
        return tmp
    else:
        return ''

    #TODO: remove below once above function is confirmed.
    # this is unnecessary and fails for very short common names
    # # find longest common substring
    # # based on https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring
    # m = [[0] * (1 + len(file2)) for i in range(1 + len(file1))]
    # longest, x_longest = 0, 0
    # for x in range(1, 1 + len(file1)):
    #     for y in range(1, 1 + len(file2)):
    #         if file1[x - 1] == file2[y - 1]:
    #             m[x][y] = m[x - 1][y - 1] + 1
    #             if m[x][y] > longest:
    #                 longest = m[x][y]
    #                 x_longest = x
    #         else:
    #             m[x][y] = 0
    # return file1[x_longest - longest: x_longest]


def find_base(read_file, read_files):
    # split the sample name from illumina tags
    replicate = '_'.join(read_file.split('_')[:-4])

    longest = 0
    rep = ''
    for read in read_files:
        tmp = lcs(replicate, read)
        if len(tmp) > longest:
            longest = len(tmp)
            rep = tmp
    return rep.strip('_')

def print_experiment(groups):
    pstr = ''
    for group, replicates in groups.items():
        pstr += group + '\n'
        for replicate, reads in replicates.items():
            pstr += '\t' + replicate + '\n'
            reads.sort()
            for read in reads:
                pstr += '\t\t' + read + '\n'
    return pstr

def find_groups(reads):

    # check if there are reads with the same name under different folders
    if not check_read_uniqueness(reads):
        return False

    data_files = {read: dirName for dirName, read_files in reads.items() for read in read_files}

    groups = dict()
    found = set()
    for read, dirName in data_files.items():
        replicate = find_base(read, data_files.keys())
        basename = find_base(read,[data_file for data_file in data_files.keys() if not data_file.startswith(replicate)])
        found.add(basename)
        if not groups.get(basename):
            groups[basename] = {replicate: ['/'.join([dirName, read])]}
        elif groups.get(basename) and not groups[basename].get(replicate):
            groups[basename][replicate] = ['/'.join([dirName, read])]
        elif groups.get(basename) and groups[basename].get(replicate):
            groups[basename][replicate].append('/'.join([dirName, read]))

    logger.info('Experiment setup:\n{}'.format(print_experiment(groups)))
    return groups




    #     pairs = []
    #     full_path = []
    #     easy_read = []
    #
    #     for forward_file in glob.glob(os.path.join(my_input, '*_1.fastq.gz')):
    #         forward_path, forward_name = forward_file.rsplit("/", 1)
    #         sample_id, ext = forward_name.rsplit("_", 1)
    #         reverse_name = sample_id + '_2.fastq.gz'
    #
    #         if os.path.isfile(reverse_name) is True:
    #             pairs.append((forward_name, reverse_name))
    #             full_path.append((forward_file, os.path.join(forward_path+"/"+reverse_name)))
    #
    #     for index, item in enumerate(pairs, start=1):
    #         easy_read.append((index, item))
    #
    #     table = tabulate(easy_read, headers=["#", "Pair"], tablefmt="grid")
    #     print (table)
    #
    #     mates_q = input("Are these pairings correct? ([y]/n) ")
    #     if mates_q.lower() != "y":
    #         print("Run aborted. Double check file names.")
    #         return
    #     else:
    #         print ("Ok saving table to Matchmaker folder")
    #         filename = os.path.join(output_path, 'table_pairs.txt')
    #         f = open(filename, 'w')
    #         f.write(table)
    #         f.close()
    #
    #     # correct delimiters for the cluster
    #     data_files = []
    #     for i in range(len(full_path)):
    #         data_files.append(str.join(' ', full_path[i]))
    #
    # else:
    #     print("Run aborted.")
    #     return
    #
    # task_count = len(data_files)
    # assert task_count > 0, "Could not find any fastq files in folder %s" % my_input
    #
    #
    # list_index=os.path.join(output_path,'read_files.txt')
    # f=open(list_index,'w')
    # f.write(str(data_files))
    # f.close()
    # print ("All done! Files are located in %s" %output_path)

# get path and run!

