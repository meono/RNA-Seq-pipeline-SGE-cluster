#!/usr/bin/env python
# vim: set fileencoding=utf-8
"""Module docstring
Template version: 1.2
"""

# TO DO: enable run threads

# for python2
from __future__ import division, print_function
import argparse
import sys
import os
import logging
import glob
import tabulate
from tabulate import tabulate

VERSION = '%(prog)s 1.01'

# for interactive call: do not add multiple times the handler
if 'LOG' not in locals():
    LOG = None
LOG_LEVEL = logging.ERROR
FORMATER_STRING = ('%(asctime)s - %(filename)s:%(lineno)d - '
                   '%(levelname)s - %(message)s')


def create_path_if_not_exists(path): # self-explanatory
    if not os.path.exists(path):
        os.makedirs(path)


def configure_log(level=LOG_LEVEL, log_file=None):
    """Configure logger
    :param log_file:
    :param level:
    """
    if LOG:
        LOG.setLevel(level)
        return LOG
    log = logging.getLogger('%s log' % os.path.basename(__file__))
    if log_file:
        handler = logging.FileHandler(filename=log_file)
    else:
        handler = logging.StreamHandler(sys.stderr)
    log_formatter = logging.Formatter(FORMATER_STRING)
    handler.setFormatter(log_formatter)
    log.addHandler(handler)
    log.setLevel(level)
    return log
LOG = configure_log()


def create_parser():
    """Return the argument parser"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--name', dest='name', required=True,
                        help='The name of this project.')

    parser.add_argument('-b', '--base', dest='base_path', default='',
                        help='''Optional base path that will be used for all other paths. This should be absolute path
                        Useful if all other paths are subfolders for the base path.''')

    parser.add_argument('-i', '--input', dest='input_path', required=True,
                        help='''The input folder path, containing the data,  absolute path if no base path was provided.
                         All files matching the the *.gz pattern will be used.''')

    parser.add_argument('-o', '--output', dest='output_path', required=True,
                        help='The output folder, where outcomes of each step will go in absolute path')

    parser.add_argument('-s', '--step', dest='step', required=True)

    # Generic arguments
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help=argparse.SUPPRESS)
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument('-q', '--quiet', '--silent', dest='quiet',
                           action='store_true', default=False,
                           help='run as quiet mode')
    verbosity.add_argument('-v', '--verbose', dest='verbose',
                           action='store_true', default=False,
                           help='run as verbose mode')
    return parser


def invoke_cluster(path): # submits job to cluster (qsub command) and gives permission for the folders created
    print("Running jobs on cluster with command 'qsub %s'. FOR SCIENCE!" % path)
    os.chmod(path, 0o755)
    os.system("qsub %s" % path)
    print("Success.")


def write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step):
    # writes a bash script as needed by the Sun Grid Engine cluster
    bash_script = '''

#!/bin/bash
#
#$ -S /bin/bash
#$ -o %(output_path)s
#$ -e %(output_path)s
#$ -r y
#$ -j y
#$ -l mem_free=%(mem_req)s
#$ -l arch=linux-x64
#$ -l netapp=2G,scratch=5G
#$ -l h_rt=%(time_req)s
#$ -t 1-%(task_count)s


inputs=(0 %(data_joined)s)
input=${inputs[$SGE_TASK_ID]}
OUT="%(output_path)s"

echo '======================================================================================================='
echo "The number of slots is:" $NSLOTS
echo "Job ID is:" $JOB_ID
echo "SGE Task ID:" $SGE_TASK_ID
echo '======================================================================================================='
echo "Input for this task: " $input
echo "Output goes to:" %(output_path)s
echo "You are at step:" %(step)s
echo '======================================================================================================='

hostname
date

%(command)s

date
qstat -j $JOB_ID

    ''' % {'output_path': output_path, 'task_count': task_count, 'data_joined': str.join(' ', data_files),
                  'mem_req': mem_req, 'time_req': time_req, 'step': step, 'command': command}

    print('========================================================================================================\n')
    print(bash_script)
    print('\n========================================================================================================')
    answer = raw_input("How you feel 'bout submitting that to the cluster? [y/n] ")

    if answer.lower() != 'y':
        print('Nevermind!')
        return

    bash_path = os.path.join(output_path, name+'_'+'script.sh')

    text_file = open(bash_path, "w")
    text_file.write(bash_script)
    text_file.close()

    invoke_cluster(bash_path)


def run_fastqc(name, input_path, output_path, step):
    # runs a quality check on fastq files

    qc_check=raw_input("Are you QCing fastq or BAM files? [FASTQ/BAM]")

    if qc_check.lower() == "fastq":
        data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))

    elif qc_check.lower() == "bam":
        data_files = glob.glob(os.path.join(input_path, '*.bam'))

    else:
        print("Not an accepted file format. Run Aborted")
        return

    task_count = len(data_files)
    print (task_count)
    assert task_count > 0, "Could not find any FASTQ .gz files in folder %s" % input_path

    output_path = os.path.join(output_path, '1.FASTQC')
    create_path_if_not_exists(output_path)

    mem_req = "10G"
    time_req = "05:00:00"
    command = 'fastqc $input --outdir=$OUT'

    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)


def run_fastx_trimmer(name, input_path, output_path, step):
    # trims 10 first bases of each read, returns files with "trimmed" prefix in chosen folder

    data_files = glob.glob(os.path.join(input_path, '*.gz'))
    task_count = len(data_files)
    assert task_count > 0, "Could not find any .gz files in folder %s" % input_path

    output_path = os.path.join(output_path, '2.FASTX')
    create_path_if_not_exists(output_path)

    mem_req = "10G"
    time_req = "10:00:00"

    command = "gzip -cd $input | fastx_trimmer -f10 -Q33 | gzip -c > $OUT/trimmed_${input##*/}"

    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)


def run_star(name, input_path, output_path, step):
    paired_q = raw_input("Do you have paired reads in separate files? [y/n] ")

    if paired_q.lower() == 'n':
        print('Ok, one file per task.')

        data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))

    elif paired_q.lower() == 'y':
        print ('Ok, looking for read mates. Mates must be labeled *_1.fastq.gz and *_2.fastq.gz for pairing.')

        # build forward and reverse file lists, match pairs
        forward = []
        reverse = []
        for read_file in glob.glob(os.path.join(input_path, '*_2.fastq.gz')):
            reverse.append(read_file)
        for read_file in glob.glob(os.path.join(input_path, '*_1.fastq.gz')):
            forward.append(read_file)

        pairs = []
        for forward_file in forward:
            for reverse_file in reverse:
                for_path, ext_for = forward_file.rsplit("_", 1)
                rev_path, ext_rev = reverse_file.rsplit("_", 1)
                if for_path == rev_path:
                    pairs.append((forward_file, reverse_file))

        # double check correct pairing in  easy to read table
        easy_read = []
        for i in range(len(pairs)):
            path1, name1 = pairs[i][0].rsplit("/", 1)
            path2, name2 = pairs[i][1].rsplit("/", 1)
            easy_read.append((i+1, name1, name2))

        table = tabulate(easy_read, headers=["Pair #", "Forward", "Reverse"], tablefmt="grid")
        print (table)
        mates_q = raw_input("Are these pairings correct? [y/n] ")
        if mates_q.lower() != "y":
            print("Run aborted.")
            return
        else:
            filename = os.path.join(output_path, 'table_pairs.txt')
            f = open(filename, 'w')
            f.write(table)
            f.close()

        data_paired = []
        for i in range(len(pairs)):
            data_paired.append(str.join(' ', pairs[i]))
        data_files = data_paired

    else:
        print("Run aborted.")
        return

    task_count = len(data_files)
    assert task_count > 0, "Could not find any fastq files in folder %s" % input_path

    command = 'input1=${inputs[$SGE_TASK_ID+$SGE_TASK_ID-1]}\n' \
    'input2=${inputs[$SGE_TASK_ID+$SGE_TASK_ID]} \n' \
    'echo "Actual input for this task is:" $input1 $input2 \n' \
    'STAR -runThread %somethings --genomeDir /netapp/home/dreuxj/Annotation/GRCh38_Gencode24/ --readFilesIn $input1 $input2\
    --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonical\
    --outFilterType BySJout --outFileNamePrefix $OUT/$SGE_TASK_ID'

    output_path = os.path.join(output_path, '3.STAR')
    create_path_if_not_exists(output_path)
    mem_req = "35G"
    time_req = "99:00:00"
    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)


def run_samtools(name, input_path, output_path, step):
    # samtools is installed on the cluster for all users but it's an older version
    # - need to ref whole path for my more recent version

    output_path = os.path.join(output_path, '4.SAMTOOLS')
    create_path_if_not_exists(output_path)

    time_req = "48:00:00"
    mem_req = "10G"

    if step == "samtools_view":
            
        data_files = glob.glob(os.path.join(input_path, '*.sam'))
        command = "/netapp/home/dreuxj/bin/samtools view -b $input -o $OUT/aligned.bam "

    elif step == "samtools_sort":

        data_files = glob.glob(os.path.join(input_path, '*.bam'))
        command = "/netapp/home/dreuxj/bin/samtools sort -o $OUT/sorted.aligned.bam $input"

    elif step == "samtools_idx":

        data_files = glob.glob(os.path.join(input_path, '*.bam'))
        command = "/netapp/home/dreuxj/bin/samtools index $input"

    elif step == "samtools_stats":

        data_files = glob.glob(os.path.join(input_path, '*.bam'))
        command = "/netapp/home/dreuxj/bin/samtools flagstat $input > ${input}_flagstats;\
         /netapp/home/dreuxj/bin/samtools idxstats $input > ${input}_idxstats"

    else:
        return

    task_count = len(data_files)
    assert task_count > 0, "Could not find any bam files in folder %s" % input_path
    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)


def run_htseq(name, input_path, output_path, step):
    # here think carefully about strand settings - use IGV if needed (not on all reads)

    data_files = glob.glob(os.path.join(input_path, '*.sorted.aligned.bam'))
    task_count = len(data_files)
    assert task_count > 0, "Could not find any sorted bam files in folder %s" % input_path

    command = 'python -m HTSeq.scripts.count -f bam -r pos -i gene_name -q -s no $input\
     /netapp/home/dreuxj/Annotation/GRCh38_Gencode24/gencode.v24.primary_assembly.annotation.gtf'

    output_path = os.path.join(output_path, '5.HTSEQ')
    create_path_if_not_exists(output_path)

    mem_req = "5G"
    time_req = "99:00:00"
    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)


def main(argv=None):
    """Program wrapper
    :param argv:
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.quiet:
        LOG.setLevel(logging.CRITICAL)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    name = args.name
    step = args.step
    input_path = os.path.join(args.base_path, args.input_path)
    output_path = os.path.join(args.base_path, args.output_path)

    create_path_if_not_exists(output_path)

    if step == 'fastqc':
        run_fastqc(name, input_path, output_path, step)
    elif step == 'fastx':
        run_fastx_trimmer(name, input_path, output_path, step)
    elif step == 'star':
        run_star(name, input_path, output_path, step)
    elif step == 'samtools_view' or step == 'samtools_sort' or step == 'samtools_idx' or step == 'samtools_stats':
        run_samtools(name, input_path, output_path, step)
    elif step == 'htseq':
        run_htseq(name, input_path, output_path, step)

    else:
        LOG.error('Did not understand step "%s". Possible values are as follows:\
         fastqc, fastx, star, htseq,samtools_view, samtools_sort, samtools_idx, samtools_stats, Run aborted.' % step)
        return 1
    return 0

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    sys.exit(main())
