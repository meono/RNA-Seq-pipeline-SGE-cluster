from __future__ import division, print_function
import os
import logging
import iLoop_RNAseq_pipeline
import readline
import pandas as pd
from validate_email import validate_email

logger = logging.getLogger(__name__)

# get tab-completion for inputs
readline.parse_and_bind("tab: complete")
readline.set_completer_delims(' \t\n;')


def get_defaults():
    """Get default parameters and setting.
    These setting are pre-defined for common analysis. For an alternate setup, generate \"RNAseq_pipeline_defaults.txt\" under home folder and/or run path.
    Pre-defined parameters are superseded by \"RNAseq_pipeline_defaults.txt\" file under home path and that, in turn, is superseded by the one under run path.
    """

    # get package defaults
    with open(os.path.join(iLoop_RNAseq_pipeline.__path__[0], 'defaults', 'RNAseq_pipeline_defaults.txt')) as rpd:
        defaults = {}
        for line in rpd.readlines():
            if line.strip():
                defaults[line.split(',')[0].strip()] = line.split(',')[1].strip()

    try:
        with open(os.path.join(os.path.expanduser("~"), 'RNAseq_pipeline_defaults.txt')) as rpd:
            for line in rpd.readlines():
                if line.strip():
                    defaults[line.split(',')[0].strip()] = line.split(',')[1].strip()
    except FileNotFoundError:
        logger.warning('"RNAseq_pipeline_defaults.txt" does not exist under home path. An email address and project ID should be should be define in that file.')

    # replace with user defaults
    try:
        with open('RNAseq_pipeline_defaults.txt') as rpd:
            for line in rpd.readlines():
                if line.strip():
                    defaults[line.split(',')[0].strip()] = line.split(',')[1].strip()
    except FileNotFoundError:
        logger.info(
            '"RNAseq_pipeline_defaults.txt" does not exist under this folder. Defaults from the package and home path will be used.')

    if 'email' not in defaults:
        if not validate_email(defaults['email']):
            while True:
                email = input('Enter a valid email address for job status: \n')
                if validate_email(email):
                    defaults['email'] = email
                    print('Writing email to "RNAseq_pipeline_defaults.txt" under home path.')
                    f = open(os.path.join(os.path.expanduser("~"), 'RNAseq_pipeline_defaults.txt'), 'w+')
                    f.write('\nemail,{}'.format(email))
                    f.close()
                    break
                else:
                    print('{} is not valid, try again.'.format(email))

    if ('project' not in defaults) or (defaults['project'] == 'projectid'):
        project = input('Enter Computerome project ID for billing: \n')
        # TODO It is possible to validate this by checking folder name under "/home/projects".
        defaults['project'] = project
        print('Writing project ID to "RNAseq_pipeline_defaults.txt" under home path.')
        f = open(os.path.join(os.path.expanduser("~"), 'RNAseq_pipeline_defaults.txt'), 'w+')
        f.write('\nproject,{}'.format(project))
        f.close()

    while True:
        # TODO: ideally this should check whether the given virtual environment has the package installed.
        # try:
        #     import iLoop_RNAseq_pipeline
        #     break
        # except ImportError:
        #     logger.warning('iLoop_RNAseq_pipeline is not available.')
        #     pass

        if (('RNAseq_env' not in defaults) or (defaults['RNAseq_env'] == 'python_environment_for_RNAseq')):
            if os.getenv('RNAseq_env') is not None:
                defaults['RNAseq_env'] = os.getenv('RNAseq_env')
                break
            else:
                defaults['RNAseq_env'] = input('Enter python environment for RNAseq where iLoop_RNAseq_pipeline is installed.')
                break

    return defaults


def get_reference(strain_code, project_path):
    ref_fs = [os.getcwd(), project_path, os.path.expanduser('~'),
              os.path.join(iLoop_RNAseq_pipeline.__path__[0], 'defaults')]
    strain_codes={}
    for ref_f in ref_fs:
        if os.path.isfile(os.path.join(ref_f, 'RNAseq_pipeline_references.tsv')):
            refs = pd.read_table(os.path.join(ref_f, 'RNAseq_pipeline_references.tsv'))
            if strain_code in list(refs.code.unique()):
                ref = refs[refs.code == strain_code].to_dict(orient='records')[0]
                logger.info('"RNAseq_pipeline_references.tsv" file found under {}'.format(ref_f))
                return ref
                strain_codes = strain_codes.union(set(refs.code.unique()))

    logger.info('"RNAseq_pipeline_references.tsv" does not contain strain code {}. Strain codes currently available: \n {}'.format(strain_code, ', '.join(list(strain_codes))))
    return False


def check_read_paths(read_path):
    if isinstance(read_path, str):
        read_path = set(path.strip() for path in read_path.split(','))
    elif isinstance(read_path, list):
        read_path = set(read_path)
    else:
        logger.warning('Read path should be a comma separated strings or list of strings for paths containing reads.')
        return False

    paths2remove = set()
    for path in read_path:
        if not os.path.exists(path):
            logger.info('{} does not exist, removed from list of paths'.format(path))
            paths2remove.add(path)

    read_path = read_path.difference(paths2remove)
    if not read_path:
        logger.warning('Read path(s) does not exist.')
        return False

    reads = dict()
    for path in read_path:
        for dirName, subdirList, fileList in os.walk(path):
            for file in fileList:
                if file.endswith('.fastq.gz'):
                    try:
                        reads[dirName].add(file)
                    except:
                        reads[dirName] = {file}
            if reads.get(dirName):
                logger.info('{} contains\n{}'.format(dirName, '\n'.join(['\t'+read for read in reads[dirName]])))
    if bool(reads):
        return reads
    else:
        logger.warning('No reads found under read path(s)')
        return False

def check_project_path(project_path):

    while True:
        if project_path == None:
            project_path = input('Enter path for project files or press enter to use current folder: \n')
            if project_path == '':
                project_path = os.getcwd()
            project_path = os.path.abspath(project_path)

        if not os.path.exists(project_path):
            create_path = input('Given project path does not exist. Do you want to create it? ([y]/n) ')
            if (create_path.lower() == 'y') or (create_path.lower() == ''):
                try:
                    os.mkdir(project_path)
                    return project_path
                except PermissionError:
                    print('Can not create folder there: Permission error. Try again')
                    project_path = input('Enter new path for project output: \n')
                    project_path = os.path.abspath(project_path)
                except:
                    print('Can not create folder there. Try again')
                    project_path = input('Enter new path for project output: \n')
                    project_path = os.path.abspath(project_path)

        if os.path.exists(project_path) and (os.access(project_path, os.W_OK)):
            if (os.listdir(project_path)):
                logger.warning('Files under project path may be overwriten.')
            return project_path
        elif os.path.exists(project_path) and not (os.access(project_path, os.W_OK)):
            logger.warning('You do not have write permission in this folder.')
            project_path = input('Enter new path for project output: \n')
            project_path = os.path.abspath(project_path)

def set_project(project_path=None, read_path=None):

    project_path = check_project_path(project_path)
    os.chdir(project_path)

    if read_path:
        reads = check_read_paths(read_path)
        if reads:
            return reads, project_path
        else:
            print('Read path(s) does not contain reads.')
            read_path = None

    if read_path == None:
        reads = check_read_paths(project_path)
        if reads:
            use_reads = input('Reads found under project folder. Do you want to use them? ([y]/n) ')
            if (use_reads.lower() == 'y') or (use_reads.lower() == ''):
                multiple_paths = input('Are there other reads to use? ([y]/n) ')
            elif use_reads.lower() == 'n':
                return reads, project_path
            else:
                logger.error('That was not y or n. Start over')
                return False

    while True:
        try:
            if (multiple_paths.lower() == 'y') or (multiple_paths.lower() == ''):
                new_path = input('Enter path for new reads: (press enter to stop)\n')
                new_path = os.path.abspath(new_path)
                if new_path == '':
                    break
                if new_path in reads.key():
                    logger.info('This folder is already added.')
                    continue
                new_reads = check_read_paths(new_path)
                if new_reads != False:
                    for newpath, newread in new_reads.items():
                        reads[newpath] = newread
                elif new_reads == False:
                    logger.warning('No reads found there.')
                    continue
            elif multiple_paths.lower() == 'n':
                break
            else:
                logger.error('That was not y or n. Start over')
                return False
        except:
            pass

        read_path = input('Enter path for reads:\n')
        reads = check_read_paths(os.path.abspath(read_path))

        if reads:
            multiple_paths = input('Are there other reads to use? ([y]/n) ')
            continue
        elif reads == False:
            logger.warning('No reads found there.')
            continue

    if bool(reads):
        return reads, project_path
    else:
        return False