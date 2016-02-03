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

VERSION = '%(prog)s 1.0'

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


myLOG = configure_log()


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
    # subprocess.Popen("qsub %s"%path)
    os.system("qsub %s" % path)
    print("Success. Cancer solved.")


def run_fastqc(name, input_path, output_path, tools_path):
    print('Setting up fastqc analysis...')

    # Check that everything is in place

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

    do_this_job = make_a_script(name, input_path, output_path, tools_path)


def make_a_script(name, input_path, output_path, tools_path):
    bash_script = '''#!/bin/bash
#                                  
#$ -S /bin/bash                    
#$ -o %(output_path)s                         #-- output directory (fill in)
#$ -e %(output_path)s                         #-- error directory (fill in)
#$ -r y                                 #-- tell the system that if a job crashes, it should be restarted
#$ -j y                                 #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=0.5G                       #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64                    #-- SGE resources (CPU type)
#$ -l netapp=1G,scratch=1G              #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                     #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-%(task_count)o                  #-- number of tasks if desired (see Tips section)

# Anything under here can be a bash script

# If you used the -t option above, this same script will be run for each task,
# but with $SGE_TASK_ID set to a different value each time (1-10 in this case).
# The commands below are one way to select a different input (PDB codes in
# this example) for each task.  Note that the bash arrays are indexed from 0,
# while task IDs start at 1, so the first entry in the tasks array variable
# is simply a placeholder

inputs=(0 %(data_joined)s)
input=${inputs[$SGE_TASK_ID]}

echo "Job ID is:" $JOB_ID
echo "SGE Task ID:" $SGE_TASK_ID
echo "Input for this task: " $input

hostname
date
%(fastqc_path)s $input --outdir=%(output_path)s
date

qstat -j $JOB_ID                                  
    ''' % {'output_path': output_path, 'task_count': len(data_files), 'data_joined': str.join(' ', data_files),
           'fastqc_path': fastqc_path}

    print('========================================================================================================\n')
    print(bash_script)
    print('\n========================================================================================================')
    answer = raw_input("How you feel 'bout submitting that to the cluster? [y/n] ")

    if answer.lower() != 'y':
        print('K, bye!')
        return

    bash_path = os.path.join(output_path, name + 'fastqc_submit_script.sh')

    text_file = open(bash_path, "w")
    text_file.write(bash_script)
    text_file.close()

    invoke_cluster(bash_path)


def main(argv=None):
    """Program wrapper
    :param argv:
    """
    if argv is None:
        argv = sys.argv[1:]
    parser = create_parser()
    args = parser.parse_args(argv)

    # Set Logging levels
    # Logging basics: 
    # Use as much as possible to see what's going on. Run with --debug, -q or -v flags to see effects.
    # LOG.debug('Debug log.')
    # LOG.info('This is a verbose log')
    # LOG.warn('This is worrisome')
    # LOG.error('This is bad')
    # LOG.fatal('Everything is ruined')
    # print("Normal output")\

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
    elif step == 'fastx_quality_stats':
        run_fastx_quality_stats(args.name, input_path, output_path, tools_path)
    else:
        LOG.error('Did not understand step "%s". Possible values are fastqc, fastqc_post. Run aborted.' % step)
        return 1

    return 0


if __name__ == '__main__':
    import doctest

    doctest.testmod()
    sys.exit(main())
