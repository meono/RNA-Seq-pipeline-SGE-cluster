from __future__ import division, print_function
import collections
import json
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
    # find longest common substring from *beginning of the string*
    for x in range(1, 1 + len(file1)):
        if file1[0:x] == file2[0:x]:
            tmp = file1[0:x]
    if 'tmp' in locals():
        return tmp
    else:
        return ''


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
    logger.info('Writing groups to json: groups.json')
    json.dump(groups, open('groups.json', 'w'))
    return groups

