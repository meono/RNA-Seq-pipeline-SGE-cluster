__author__ = 'emre'

import os
import stat
from setuptools import setup

# run_RNAseq file needs to be executable but git won't set that properly
# os.chmod(os.path.join(os.getcwd(), 'iLoop_RNAseq_pipeline', 'run_RNAseq.py'),
#          os.stat(os.path.join(os.getcwd(), 'iLoop_RNAseq_pipeline',
#                               'run_RNAseq.py')).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

requirements = ['validate_email>=1.3',
                'pandas>=0.18.1']

setup(
    name='iLoop_RNAseq_pipeline',
    version='0.1dev',
    packages=['iLoop_RNAseq_pipeline'],
    license='Apache License, Version 2.0',
    include_package_data=True,
    install_requires=requirements,
    # data_files=[('iLoop_RNAseq_pipeline', ['iLoop_RNAseq_pipeline/defaults/References.tsv', 'iLoop_RNAseq_pipeline/defaults/RNAseq_pipeline_defaults.txt'])],
    # long_description=open('dontREADMEyet.md').read(),
)
