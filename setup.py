__author__ = 'emre'

from setuptools import setup

requirements = ['validate_email>=1.3']

setup(
    name='iLoop_RNAseq_pipeline',
    version='0.1dev',
    packages=['iLoop_RNAseq_pipeline',],
    license='Apache License, Version 2.0',
    include_package_data=True,
    install_requires=requirements,
    # data_files=[('iLoop_RNAseq_pipeline', ['iLoop_RNAseq_pipeline/defaults/References.tsv', 'iLoop_RNAseq_pipeline/defaults/RNAseq_pipeline_defaults.txt'])],
    # long_description=open('dontREADMEyet.md').read(),
)