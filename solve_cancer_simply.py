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
import functools
import logging
import glob
import operator
from functools import reduce
import subprocess

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
                        help='The name of this project (e.g. VS1291).')

    parser.add_argument('-b', '--base', dest='base_path', default='',
                        help='''Optional base path that will be used for all other paths. This should be absolute path
                        Useful if all other paths are subfolders for the base path.''')

    parser.add_argument('-i', '--input', dest='input_path', required=True,
                        help='''The input folder path, containing the data,  absolute path if no base path was provided.
                         All files matching the the *.gz pattern will be used.''')

    parser.add_argument('-o', '--output', dest='output_path', required=True,
                        help='The output folder, where outcomes of each step will go in absolute path')

    parser.add_argument('-t', '--tools', dest='tools_path', required=True,
                        help='The tools folder, containing all the tools used by the pipeline, in absolute path')

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
    print("Success. Cancer solved.")


def write_bash_script(name, data_files, output_path, mem_req, tool_path, command, which_step):
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
#$ -l netapp=1G,scratch=1G              
#$ -l h_rt=24:00:00                     
#$ -t 1-%(task_count)o                 

inputs=(0 %(data_joined)s)
input=${inputs[$SGE_TASK_ID]}
TOOL="%(tool_path)s"
OUT="%(output_path)s"
OUTFILE=${input}_trimmed.fastq.gz
GTF_ANNOT=/netapp/home/dreuxj/hg38/Annotation/genes.gtf
BWT2_IDX=/netapp/home/dreuxj/hg38/Sequence/Bowtie2Index/genome
GENOME_FASTA=/netapp/home/dreuxj/hg38/Sequence/Bowtie2Index/genome.fa
GENOME_DIR=/netapp/home/dreuxj/GRCh38_Gencode24/

echo "Job ID is:" $JOB_ID
echo "SGE Task ID:" $SGE_TASK_ID
echo "Input for this task: " $input
echo "Output goes to:" %(output_path)s
echo "You are using the tool:" %(tool_path)s
echo "You are at step:" %(which_step)s
echo "Your command is:" %(command)s
hostname
date
%(command)s
date

qstat -j $JOB_ID

    ''' % {'output_path': output_path, 'task_count': len(data_files), 'data_joined': str.join(' ', data_files), 'mem_req': mem_req,
           'tool_path': tool_path, 'which_step': which_step, 'command': command}

    print('========================================================================================================\n')
    print(bash_script)
    print('\n========================================================================================================')
    answer = raw_input("How you feel 'bout submitting that to the cluster? [y/n] ")

    if answer.lower() != 'y':
        print('K, bye Felicia!')
        return

    bash_path = os.path.join(output_path, name+'_'+which_step+'_'+'submit_script.sh')

    text_file = open(bash_path, "w")
    text_file.write(bash_script)
    text_file.close()

    invoke_cluster(bash_path)


def run_fastqc(name, input_path, output_path, tools_path):
    print('Setting up fastqc analysis...')

    # Fast qc found?
    fastqc_path = os.path.join(tools_path, 'FastQC/fastqc')
    assert os.path.isfile(fastqc_path), "Could not find fastqc at path %s" % fastqc_path

    # Data files found?
    data_files = glob.glob(os.path.join(input_path, '*.gz'))
    assert len(data_files) > 0, "Could not find any .gz files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'fastqc')
    create_path_if_not_exists(output_path)

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    #how much memory do you need for this job?
    mem_req = "1G"

    # what is your command
    command = "$TOOL $input --outdir=$OUT"

    #what are you calling this step in the pipeline
    which_step = 'FASTQC'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, fastqc_path, command, which_step)


def run_fastx_trimmer(name, input_path, output_path, tools_path):
   
    # Tool found?
    fastx_trimmer_path = os.path.join(tools_path, 'fastx_toolkit/fastx_trimmer')
    assert os.path.isfile(fastx_trimmer_path), "Could not find trimmer at path %s" % fastx_trimmer_path

    # Data files found?
    data_files = glob.glob(os.path.join(input_path, '*.gz'))
    assert len(data_files) > 0, "Could not find any .gz files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'fastx_trimmer')
    create_path_if_not_exists(output_path)

    #how much memory do you need for this job?
    mem_req = "1G"

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    # need to unzip files for this to run.
    # what is your command
    command = "gzip -cd $input | $TOOL -f10 -Q33 | gzip -c > $OUTFILE"

    #what are you calling this step in the pipeline
    which_step = 'FASTX_TRIMMER'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, fastx_trimmer_path, command, which_step)


def run_tophat(name, input_path, output_path, tools_path):
    
    # Tool found?
    tophat_path = os.path.join(tools_path, 'tophat2/tophat2')
    assert os.path.isfile(tophat_path), "Could not find tophat2 at path %s" % tophat_path

    # Trimmed data files found? 
    data_files = glob.glob(os.path.join(input_path, '*_trimmed.fastq.gz'))
    assert len(data_files) > 0, "Could not find any -trimmed files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'tophat')
    create_path_if_not_exists(output_path)

    #how much memory do you need for this job?
    mem_req = "5G"

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    # what is your command
    # !!! run this on a SINGLE FILE FIRST. Then you can do the rest (transcriptome-idx problem)
    command = "$TOOL -G $GTF_ANNOT --transcriptome-index $OUT/trans_idx -o $OUT $BWT2_IDX $input"

    #what are you calling this step in the pipeline
    which_step = 'TOPHAT'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, tophat_path, command, which_step)

def run_STAR(name, input_path, output_path, tools_path):
    
    # Tool found?
    star_path = os.path.join(tools_path, 'STAR/STAR')
    assert os.path.isfile(star_path), "Could not find STAR at path %s" % star_path

    # Trimmed data files found? 
    data_files = glob.glob(os.path.join(input_path, '*_trimmed.fastq.gz'))
    assert len(data_files) > 0, "Could not find any -trimmed files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'STAR/')
    create_path_if_not_exists(output_path)

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    #how much memory do you need for this job?
    mem_req = "30G"

    # what is your command
    #this aligner requires 30G free RAM for the human genome. Request acordingly!
    command = "$TOOL --runThreadN 12 --genomeDir $GENOME_DIR --readFilesIn $input --readFilesCommand gunzip -c --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --outFileNamePrefix ~/JD1291/STAR/"

    #what are you calling this step in the pipeline
    which_step = 'STAR'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, star_path, command, which_step)

def run_samtools(name, input_path, output_path, tools_path):
    
    # Tool found?
    samtools_path = os.path.join(tools_path, 'samtools/Samtools')
    assert os.path.isfile(samtools_path), "Could not find samtools at path %s" % samtools_path

    data_files = glob.glob(os.path.join(input_path, '*.sam'))
    assert len(data_files) > 0, "Could not find any sam files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'Samtools')
    create_path_if_not_exists(output_path)

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    #how much memory do you need for this job?
    mem_req = "2G"

    # what is your command
    command = "$TOOL view -bS $input -o $OUT/Aligned.bam "

    #what are you calling this step in the pipeline
    which_step = 'Samtools sam to bam'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, samtools_path, command, which_step)


def run_cufflinks(name, input_path, output_path, tools_path):
    
    # Tool found?
    cufflinks_path = os.path.join(tools_path, 'cufflinks/cufflinks')
    assert os.path.isfile(cufflinks_path), "Could not find cufflinks at path %s" % cufflinks_path    

    data_files = glob.glob(os.path.join(input_path, '*.sorted.bam'))
    assert len(data_files) > 0, "Could not find any sorted bam files in folder %s" % input_path

    # Setup the output for this step
    output_path = os.path.join(output_path, 'Cufflinks')
    create_path_if_not_exists(output_path)

    # Compute size of input (can be useful for runtime limits, below).
    sizes = list(map(os.path.getsize, data_files))
    total_size = reduce(operator.add, sizes)

    #how much memory do you need for this job?
    mem_req = "15G"

    # what is your command
    command = "$TOOL -m 42 -vu -G $GTF_ANNOT -b $GENOME_FASTA $input -o $OUT"

    #what are you calling this step in the pipeline
    which_step = 'Cufflinks'

    # send to bash script function
    write_bash_script(name, data_files, output_path, mem_req, cufflinks_path, command, which_step)


def run_cuffdiff(name, input_path, output_path, tools_path):
    pass


def run_cuffmerge(name, input_path, output_path, tools_path):
    pass


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
    tools_path = os.path.join(args.base_path, args.tools_path)

    create_path_if_not_exists(output_path)

    if step == 'fastqc':
        run_fastqc(name, input_path, output_path, tools_path)
    elif step == 'trimmer':
        run_fastx_trimmer(name, input_path, output_path, tools_path)
    elif step == 'tophat':
        run_tophat(name, input_path, output_path, tools_path)
    elif step == 'STAR':
        run_STAR(name, input_path, output_path, tools_path)
    elif step == 'samtools':
        run_samtools(name, input_path, output_path, tools_path)    
    elif step == 'cufflinks':
        run_cufflinks(name, input_path, output_path, tools_path)
    elif step == 'cuffdiff':
        run_cuffdiff(name, input_path, output_path, tools_path)
    elif step == 'cuffmerge':
        run_cuffmerge(name, input_path, output_path, tools_path)
    else:
        LOG.error('Did not understand step "%s". Possible values are fastqc, trimmer, tophat, STAR, cufflinks, cuffdiff, and\
                  cuffmerge. Run aborted.' % step)
        return 1

    return 0


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    sys.exit(main())
