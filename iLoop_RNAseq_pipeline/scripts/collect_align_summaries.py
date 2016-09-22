#!/usr/bin/env python

import pandas as pd
import json
import os
import argparse
import logging

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='This script will collect align summaries based on hisat2 run.')

parser.add_argument('-p', '--project-path', help='Path to project folder')
parser.add_argument('-g', '--groups-json', help='json file for experiment setup with full path.')
parser.add_argument('-o', '--output', help='Path and prefix of output file. e.g.: ./results/align_summaries.tsv')
parser.add_argument('-m', '--map-to-mask', help='Collect align summaries for mask only run.', action='store_true', default=False)
args = parser.parse_args()

groups = json.load(open(args.groups_json, 'r'))

if args.map_to_mask:
    pp = os.path.abspath(os.path.join(args.project_path, 'map_to_mask')) # pp: path prefix
else:
    pp = os.path.abspath(args.project_path) # pp: path prefix
for sample in [sample for group in groups.values() for sample in group.keys()]:
    with open(os.path.join(pp, sample, 'align_summary.txt'), 'r') as f:
        for line in f.readlines():
            if line.startswith('Left'):  # it is tophat summary
                sum_type = 'tophat2'
                if 'df' not in locals():
                    df = pd.DataFrame(columns=['Left_reads', 'Mapped_Left_reads', 'Percent_mapped_Left_reads',
                                               'Multiple_alignments_Left_reads',
                                               'Percent_multiple_alignments_Left_reads',
                                               'Right_reads', 'Mapped_Right_reads', 'Percent_mapped_Right_reads',
                                               'Multiple_alignments_Right_reads',
                                               'Percent_multiple_alignments_Right_reads',
                                               'Aligned_pairs', 'Multiple_alignments_pairs',
                                               'Percent_multiple_alignments_pairs', 'Concordant_pair_alignment_rate'])
            elif line.endswith('reads; of these:\n'):  # it is hisat2 summary
                sum_type = 'hisat2'
                if 'df' not in locals():
                    df = pd.DataFrame(columns=['Reads', 'Paired_reads', 'Pairs_not_concordantly_aligned',
                                               'Pairs_concordantly_aligned_once', 'Pairs_concordantly_aligned_multiple',
                                               'Pairs_discordantly_aligned_once',
                                               'Remaining_singles_not_aligned', 'Remaining_singles_aligned_once',
                                               'Remaining_singles_aligned_multiple', 'Overall_alignment_percentage'])

            if sum_type == 'tophat2':
                if line.startswith('Left'):
                    side = 'Left'
                elif line.startswith('Right'):
                    side = 'Right'
                elif line.startswith('Aligned'):
                    side = False

                if line.startswith('          Input     :') and side != False:
                    df.loc[i, side + '_reads'] = line.split()[2]
                elif line.startswith('           Mapped   :') and side != False:
                    df.loc[i, 'Mapped_' + side + '_reads'] = line.split()[2]
                    df.loc[i, 'Percent_mapped_' + side + '_reads'] = line.split('(')[1].split('%')[0]
                elif line.startswith('            of these:') and side != False:
                    df.loc[i, 'Multiple_alignments_' + side + '_reads'] = line.split()[2]
                    df.loc[i, 'Percent_multiple_alignments_' + side + '_reads'] = line.split('(')[1].split('%')[0]
                elif side == False and line.startswith('Aligned pairs:'):
                    df.loc[i, 'Aligned_pairs'] = line.split()[2]
                elif side == False and line.startswith('     of these:'):
                    df.loc[i, 'Multiple_alignments_pairs'] = line.split()[2]
                    df.loc[i, 'Percent_multiple_alignments_pairs'] = line.split('(')[1].split('%')[0]
                elif side == False and line.endswith('concordant pair alignment rate.\n'):
                    df.loc[i, 'Concordant_pair_alignment_rate'] = line.split('%')[0]

            elif sum_type == 'hisat2':
                if line.endswith('reads; of these:\n'):
                    df.loc[i, 'Reads'] = line.split()[0]
                elif line.endswith('were paired; of these:\n'):
                    df.loc[i, 'Paired_reads'] = line.split()[0]
                elif line.endswith('aligned concordantly 0 times\n'):
                    df.loc[i, 'Pairs_not_concordantly_aligned'] = line.split()[0]
                elif line.endswith('aligned concordantly exactly 1 time\n'):
                    df.loc[i, 'Pairs_concordantly_aligned_once'] = line.split()[0]
                elif line.endswith('aligned concordantly >1 times\n'):
                    df.loc[i, 'Pairs_concordantly_aligned_multiple'] = line.split()[0]
                elif line.endswith('aligned discordantly 1 time\n'):
                    df.loc[i, 'Pairs_discordantly_aligned_once'] = line.split()[0]
                elif line.endswith('aligned 0 times\n'):
                    df.loc[i, 'Remaining_singles_not_aligned'] = line.split()[0]
                elif line.endswith('aligned exactly 1 time\n'):
                    df.loc[i, 'Remaining_singles_aligned_once'] = line.split()[0]
                elif line.endswith('aligned >1 times\n'):
                    df.loc[i, 'Remaining_singles_aligned_multiple'] = line.split()[0]
                elif line.endswith('overall alignment rate\n'):
                    df.loc[i, 'Overall_alignment_percentage'] = line.split('%')[0]


df.to_csv('align_summaries.csv', sep='	')
