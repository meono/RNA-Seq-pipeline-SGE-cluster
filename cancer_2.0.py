#!/usr/bin/env python
# vim: set fileencoding=utf-8
"""Module docstring
Template version: 1.2
"""

# for python2
from __future__ import division, print_function

import argparse
import sys
import os
import logging
import glob

VERSION = '%(prog)s 1.01'

# for interactive call: do not add multiple times the handler
if 'LOG' not in locals():
    LOG = None
LOG_LEVEL = logging.ERROR
FORMATER_STRING = ('%(asctime)s - %(filename)s:%(lineno)d - '
                   '%(levelname)s - %(message)s')


def create_path_if_not_exists(path):
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


def invoke_cluster(path):
    print("Running jobs on cluster with command 'qsub %s'. FOR SCIENCE!" % path)
    os.chmod(path, 0o755)
    os.system("qsub %s" % path)
    print("Success. Cancer 2.0 solved.")


def write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step):
    bash_script = '''

#!/bin/bash
#
#$ -S /bin/bash
#$ -o {output_path:s}
#$ -e {output_path:s}
#$ -r y
#$ -j y
#$ -l mem_free={mem_req:s}
#$ -l arch=linux-x64
#$ -l netapp=2G,scratch=5G
#$ -l h_rt={time_req:o}
#$ -t 1-{task_count:o}

inputs=(0 {data_joined:s})
input=${{inputs[$SGE_TASK_ID]}}



echo '======================================================================================================='
echo "Job ID is:" $JOB_ID
echo "SGE Task ID:" $SGE_TASK_ID
echo '======================================================================================================='
echo "Input for this task: " $input
echo "Output goes to:" {output_path:s}
echo "You are at step:" {step:s}
echo "Your command is:" {command:s}
echo '======================================================================================================='

hostname
date

{command:s}

date

qstat -j $JOB_ID

    '''.format(**{'output_path': output_path, 'task_count': task_count, 'data_joined': str.join(' ', data_files),
                  'mem_req': mem_req,
                  'time_req': time_req, 'step': step, 'command': command})

    print('========================================================================================================\n')
    print(bash_script)
    print('\n========================================================================================================')
    answer = raw_input("How you feel 'bout submitting that to the cluster? [y/n] ")

    if answer.lower() != 'y':
        print('K, bye Felicia!')
        return

    bash_path = os.path.join(output_path, name+'_'+step+'_'+'_script.sh')

    text_file = open(bash_path, "w")
    text_file.write(bash_script)
    text_file.close()

    invoke_cluster(bash_path)

base_path="/netapp/home/dreuxj/"

def run_fastqc(name, input_path, output_path, step):

    data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))
    task_count=len(data_files)
    assert task_count> 0, "Could not find any FASTQ .gz files in folder %s" % input_path

    output_path = os.path.join(output_path, '1.Fastqc')
    create_path_if_not_exists(output_path)

    mem_req = "5G"
    time_req ="05:00:00"
    command = "fastqc $input --outdir={output_path:s}"

    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

def run_fastx_trimmer(name, input_path, output_path, step):

    data_files = glob.glob(os.path.join(input_path, '*.gz'))
    task_count=len(data_files)
    assert task_count > 0, "Could not find any .gz files in folder %s" % input_path

    output_path = os.path.join(output_path, '2.fastx_trimmer')
    create_path_if_not_exists(output_path)

    mem_req = "10G"
    time_req="10:00:00"
    command = "gzip -cd $input | fastx_trimmer -f10 -Q33 | gzip -c > ${input}_trimmed.fastq.gz"

    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

def run_star(name, input_path, output_path, step):

    data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))
    task_count = len(data_files)
    assert task_count > 0, "Could not find any fastq files in folder %s" % input_path

    output_path = os.path.join(output_path, '3.STAR')
    create_path_if_not_exists(output_path)


    mem_req = "30G"
    time_req="48:00:00"
    command = "STAR --runThreadN 12 --genomeDir /netapp/home/dreuxj/GRCh38_Gencode24/ --readFilesIn $input \
    --readFilesCommand gunzip -c --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated\
     --outFilterType BySJout --outFileNamePrefix {output_path:s}/${input}"

    write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

def run_samtools(name, input_path, output_path, step):

    output_path = os.path.join(output_path,step)
    create_path_if_not_exists(output_path)

    data_files = glob.glob(os.path.join(input_path, '*.bam'))
    task_count = len(data_files)
    assert task_count > 0, "Could not find any bam files in folder %s" % input_path

    time_req="20:00:00"
    mem_req="10G"

    if step == "samtools_view":
            
        data_files = glob.glob(os.path.join(input_path, '*.sam'))
        command = "samtools view -bS $input -o $OUTPUT/${input}.aligned.bam "
        write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

    elif step == "samtools_sort":

        command = "samtools sort $input -o $OUTPUT/${input}.sorted.bam"
        write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

    elif step == "samtools_idx":

        command = "samtools index $input"
        write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

    elif step == "samtools_stats":

        command = "samtools flagstat $input > ${input}_flagstats"
        "samtools idxstats $input > ${input}_idxstats"
        write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)

    else:
        print("Can't match step to function, run aborted")


def run_cufflinks_suite(name, input_path, output_path, step):

    output_path = os.path.join(output_path, step)
    create_path_if_not_exists(output_path)
    data_files=""
    task_count =len(data_files)


    if step == "cufflinks":

        data_files = glob.glob(os.path.join(input_path, '*.sorted.bam'))
        mem_req = "30G"
        command = "GTF_ANNOT=/netapp/home/dreuxj/hg38/Annotation/genes.gtf" \
                  "GENOME_FASTA=/netapp/home/dreuxj/hg38/Sequence/Bowtie2Index/genome.fa" \
                  "cufflinks -m 42 -u -G $GTF_ANNOT -b $GENOME_FASTA $input -o $OUT"

    elif step == "cuffdiff":

        data_files = glob.glob(os.path.join(input_path, '*.sorted.bam'))

        mem_req = "10G"
        command = "RECTUS_BAM=/netapp/home/dreuxj/JD1291/Samtools/Rectus/Aligned_rectus.sorted.bam" \
                "VASTUS_BAM=/netapp/home/dreuxj/JD1291/Samtools/Vastus/Aligned_vastus.sorted.bam" \
                "cuffdiff -L Rectus,Vastus -b $GENOME_FASTA  -o $OUTPUT $RECTUS_BAM $VASTUS_BAM"

    elif step == "cuffmerge":

        data_files = glob.glob(os.path.join(input_path, '*.txt'))
        mem_req = "30G"
        command = "GTF_ANNOT=/netapp/home/dreuxj/hg38/Annotation/genes.gtf" \
                  "GENOME_FASTA=/netapp/home/dreuxj/hg38/Sequence/Bowtie2Index/genome.fa" \
                  "cuffmerge -g $GTF_ANNOT -s $GENOME_FASTA  -o $OUT $input"

    else:
        print("Can't match step to function, run aborted")


    assert task_count > 0, "Could not find any sorted bam or sam files in folder %s" % input_path

    write_bash_script(name, data_files, output_path, mem_req, task_count, command, step)


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
    elif step == 'samtools_view' or step =='samtools_sort' or step=='samtools_idx' or step=='samtools_stats':
        run_samtools(name, step, input_path, output_path, step)
    elif step == 'cuffdiff'or step =='cufflinks' or step=='cuffmerge':
        run_cufflinks_suite(name, input_path, output_path, step)

    else:
        LOG.error('Did not understand step "%s". Possible values are fastqc, fastx, STAR, Samtools_*, cuff*, cuffdiff.\
         Run aborted.' % step)
        return 1

    return 0


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    sys.exit(main())