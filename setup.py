__author__ = 'emre'

from setuptools import setup

requirements = ['validate_email>=1.3',
                'pandas>=0.18.1']

setup(
    name='iLoop_RNAseq_pipeline',
    version='0.1',
    packages=['iLoop_RNAseq_pipeline'],
    license='Apache License, Version 2.0',
    include_package_data=True,
    install_requires=requirements,
    long_description=open('README.md').read(),
    # data_files=[],
)
