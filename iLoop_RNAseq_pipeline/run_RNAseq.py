#!/usr/bin/env python

__author__ = 'emre'
import iLoop_RNAseq_pipeline
import iLoop_RNAseq_pipeline.initiate_project as ip
import iLoop_RNAseq_pipeline.matchmaker as mm
import iLoop_RNAseq_pipeline.master_qsub as mq
import argparse
import logging
import pandas as pd
import os

parser = argparse.ArgumentParser(description='This script will run RNAseq pipeline. Only works on Computerome.')

parser.add_argument('-p', '--project-path', help='Full path to project folder', default=None)
parser.add_argument('-r', '--read-path',
                    help='Full path to reads folder(s). Can be suplied as comma separated string. All files under the path tree will be added.',
                    default=None)
# parser.add_argument('-s', '--strain-code', help='Code for the reference strain. Available strains and corresponding codes are on "References.tsv".', required=True)
parser.add_argument('-s', '--strain-code',
                    help='Code for the reference strain. Available strains and corresponding codes are on "References.tsv".',
                    default=None)
args = parser.parse_args()

if __name__ == "__main__":
    logging.basicConfig(filename='{}/RNAseq_pipeline.log'.format(args.project_path),
                        format='%(asctime)s %(levelname)s : %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S',
                        level=logging.DEBUG)
    logging.info('\n-----\nRNAseq pipeline started\n----')
    defaults = ip.get_defaults()
    ref_fs = [os.getcwd(), project_path, os.path.expanduser('~'), os.path.join(iLoop_RNAseq_pipeline.__path__, 'defaults')]
    for ref_f in ref_fs:
        try:
            refs = pd.read_table(os.path.join(ref_f, 'RNAseq_pipeline_references.tsv'))
            ref = refs[refs.code == args.strain_code].to_dict(orient='records')[0]
            logging.info('"RNAseq_pipeline_references.tsv" file found under {}'.format(ref_f))
            break
        except FileNotFoundError:
            pass
    reads, project_path = ip.set_project(project_path=args.project_path, read_path=args.read_path)
    groups = mm.find_groups(reads)
    jobs = mq.jobsubmitter(project_path=project_path, groups=groups,readtype='raw', ref=ref, defaults=defaults, ppn=8)
