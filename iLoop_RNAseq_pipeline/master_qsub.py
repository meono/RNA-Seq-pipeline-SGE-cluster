from __future__ import division, print_function
import os
from os.path import abspath
from os.path import join as join_path
import subprocess
import logging
import sys
import iLoop_RNAseq_pipeline

logger = logging.getLogger(__name__)

job_header = \
    '''#!/bin/sh
### Account information
#PBS -W group_list=PROJECT -A PROJECT
# -- Job dependencies --
#PBS -W depend=DEPEND
# -- Name of the job --
#PBS -N JOBNAME
# -- estimated wall clock time (execution time): hh:mm:ss --
#PBS -l walltime=WALLTIME
# -- number of processors/cores/nodes --
#PBS -l nodes=1:ppn=PPN
# -- output destination --
#PBS -o JOB_OUTPUTS/JOBNAME_${PBS_JOBID%%.*}.o
#PBS -j oe
# -- user email address --
#PBS -M EMAILADDRESS
# -- run in the current working (submission) directory --
if test X$PBS_ENVIRONMENT = XPBS_BATCH; then cd $PBS_O_WORKDIR; fi
# here follow the commands you want to execute'''


def check_fastqc(groups, output_path):
    """Checks if fastqc files exists and are not empty."""
    for readf in [readf for group in groups.values() for readfs in group.values() for readf in readfs]:
        qcf = readf.split('.')[0].split('/')[-1]
        if (not os.path.isfile(join_path(output_path, qcf + '_fastqc.zip'))) or \
                (not os.path.isfile(join_path(output_path, qcf + '_fastqc.html'))):
            return False
        elif (os.path.getsize(join_path(output_path, qcf + '_fastqc.zip')) == 0) or \
                (os.path.getsize(join_path(output_path, qcf + '_fastqc.html')) == 0):
            return False
    return True


def check_mapandlink(groups, project_path):
    """Check if the cufflinks run generated HIDATA issues."""
    # TODO: This should expand to existence of previous results and possibly other data issues.
    check = True
    for group in groups.values():
        for replicate in group.keys():
            cfs = [file for file in os.listdir(join_path(project_path, replicate)) if file.endswith('tracking')]
            for file in cfs:
                with open(join_path(project_path, replicate, file), 'r') as cf:
                    s = cf.read()
                    cf.close()
                if 'HIDATA' in s:
                    check = False
                    logger.error('{} had HIDATA issue and need to be repeated with a larger bundle size.'.format(replicate))
    return check


def fastqc_job(project_path, groups, output_path, defaults, ppn='8', walltime='02:00:00'):
    """Runs fastqc for all the reads.

    reads: dictionary from set_project function
    output_path: path for the output of fastqc
                - ideally, project_path/reads/QC_output/XYZ where XYZ denotes raw or trimmed
     """

    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'fastqc')
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', '')
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs FastQC/0.11.2''']

    jobstr += ['fastqc -t PPN -o {} {}'.format(output_path,
                                               ' '.join(
                                                   [readf for group in groups.values() for readfs in group.values() for
                                                    readf in readfs]))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def mapandlink_jobs(project_path, sample, reads, defaults, ref, jobs, ppn='8', walltime='12:00:00', map_to_mask=False):
    mljobs = ['hisat2', 'stringtie', 'cufflinks', 'htseq-count', 'featureCounts']
    if map_to_mask:
        jobs = ['hisat2_to_mask']
    elif (jobs == []) and (map_to_mask == False):
        jobs = mljobs

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', '_'.join([sample] + [job for job in jobs if job in mljobs]))
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', '')
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools samtools/1.2 bowtie2/2.2.5 cufflinks/2.2.1
# hisat2 and stringtie modules are not kept up-to-date, so use the local ones:
PATH=$PATH:/home/projects/cu_10010/programs/hisat2-2.0.4:/home/projects/cu_10010/programs/stringtie-1.2.3.Linux_x86_64
export PATH''']

    jobstr += ['export hisat2_genomic_indexes={}'.format(ref['hisat2_genomic_indexes'])]
    if (map_to_mask) and (not os.path.exists(join_path(project_path, 'map_to_mask', sample))):
        jobstr += ['mkdir -p {}'.format(abspath(join_path(project_path, 'map_to_mask', sample)))]
    elif not os.path.exists(join_path(project_path, sample)):
        jobstr += ['mkdir {}'.format(abspath(join_path(project_path, sample)))]

    R1reads = [read for read in reads if read[-15:-13] == 'R1']
    R1reads.sort()
    R2reads = [read for read in reads if read[-15:-13] == 'R2']
    R2reads.sort()

    if (ref.get('hisat2_mask_indexes')) and ('hisat2_to_mask' in jobs):
        logger.info('Using hisat2 options: {}'.format(defaults['hisat2_options']))
        jobstr += ['''echo "hisat2"
    hisat2 {} -p PPN -x {} {} {} -S {} 2>{} '''.format(defaults['hisat2_options'],
                                        (ref['hisat2_mask_indexes']),
                                        ('-1' + ','.join(R1reads)),
                                        ('-2' + ','.join(R2reads)),
                                        (abspath(join_path(project_path, 'map_to_mask', sample, 'accepted_hits.sam'))),
                                        (abspath(join_path(project_path, 'map_to_mask', sample, 'align_summary.txt'))))]
    elif (ref.get('hisat2_genomic_indexes')) and ('hisat2' in jobs):
        logger.info('Using hisat2 options: {}'.format(defaults['hisat2_options']))
        jobstr += ['''echo "hisat2"
hisat2 {} -p PPN -x {} {} {} 2>{} | \
samtools view -@ PPN -hbu - | \
samtools sort -@ PPN - {}'''.format(defaults['hisat2_options'],
                                    ref['hisat2_genomic_indexes'],
                                    ('-1' + ','.join(R1reads)),
                                    ('-2' + ','.join(R2reads)),
                                    (abspath(join_path(project_path, sample, 'align_summary.txt'))),
                                    (abspath(join_path(project_path, sample, 'accepted_hits.sorted'))))]

    if 'stringtie' in jobs:
        logger.warning('Beware: Stringtie does not allow masking for now.')
        logger.info('Using stringtie options: {}'.format(defaults['stringtie_options']))
        jobstr += ['echo "stringtie"\nstringtie {} -p PPN {} -o {} -A {} {}'.format(defaults['stringtie_options'],
                                                                                    (('-G ' + ref[
                                                                                        'gff_genome']) if ref.get(
                                                                                        'gff_genome') else ''),
                                                                                    abspath(
                                                                                        join_path(project_path,
                                                                                                     sample,
                                                                                                     'transcripts.gtf')),
                                                                                    abspath(
                                                                                        join_path(project_path,
                                                                                                     sample,
                                                                                                     'gene_abund.tab')),
                                                                                    abspath(
                                                                                        join_path(project_path,
                                                                                                     sample,
                                                                                                     'accepted_hits.sorted.bam')))]

    if 'cufflinks' in jobs:
        logger.info('Using cufflinks options: {}'.format(defaults['cufflinks_options']))
        jobstr += ['echo "cufflinks"\ncufflinks {} -p PPN {} {} -o {} {}'.format(defaults['cufflinks_options'],
                                                                                 ('-G ' + ref['gff_genome']) if ref.get(
                                                                                     'gff_genome') else '',
                                                                                 ('-M ' + ref['gff_mask']) if ref.get(
                                                                                     'gff_mask') else '',
                                                                                 (abspath(
                                                                                     join_path(project_path,
                                                                                                  sample))),
                                                                                 (abspath(
                                                                                     join_path(project_path, sample,
                                                                                                  'accepted_hits.sorted.bam'))))]

    # htseq-count doesn't work well with positional sorted bam files. Switch to featureCounts from subread package
    # if any(job for job in jobs if job in ['htseq-count', 'edgeR', 'DESeq']):
    #     logger.info('Using htseq options: {}'.format(defaults['htseq_options']))
    #     jobstr += ['echo "htseq"\nhtseq-count {} -f bam -r pos {} {} -o {} > {}'.format(defaults['htseq_options'],
    #                                                                              (abspath(
    #                                                                                  join_path(project_path, sample,
    #                                                                                               'accepted_hits.sorted.bam'))),
    #                                                                              ((ref['gff_genome']) if ref.get(
    #                                                                                  'gff_genome') else ''),
    #                                                                              (abspath(
    #                                                                                  join_path(project_path, sample,
    #                                                                                               'htseq_counts.sam'))),
    #                                                                              (abspath(
    #                                                                                  join_path(project_path, sample,
    #                                                                                               'htseq_counts.out'))))]

    # it turnsout featureCounts does a name based sorting before operating. therefore this might be just as inefficient.
    if any(job for job in jobs if job in ['featureCounts', 'edgeR', 'DESeq']):
        logger.info('Using htseq options: {}'.format(defaults['htseq_options']))
        jobstr += ['echo "featureCounts"\nfeatureCounts {} -a {} -o {} {}'.format(defaults['featureCounts_options'],
                                                                                        ((ref['gff_genome']) if ref.get(
                                                                                            'gff_genome') else ''),
                                                                                        (abspath(join_path(project_path,
                                                                                                                      sample,
                                                                                                                      'featureCounts_{}.out'.format(sample)))),
                                                                                        (abspath(join_path(project_path,
                                                                                                                      sample,
                                                                                                                      'accepted_hits.sorted.bam'))))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def collect_counts_job(project_path, output, mapjobIDs, defaults, ppn='1', walltime="00:05:00"):
    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'collect_counts')
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', (
        'afterok:{}'.format(':'.join([mapjob for mapjob in mapjobIDs])) if mapjobIDs != [''] else ''))
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    # Pass all environmental variables to the job - this should take care of the virtual environment issue
    # TODO: clear out virtual environment arguments if this works
    jobstr += ['#PBS -V']

    # this is for htseq-count, which won't be used. For now, featureCounts will be used.
    # jobstr += ['python {}/htseq_count_collector.py -p {} -g {} -o {} '.format(abspath(join_path(iLoop_RNAseq_pipeline.__path__[0], 'scripts')),
    #                                                                           abspath(project_path),
    #                                                                           abspath(join_path(join_path(project_path, 'groups.json'))),
    #                                                                           output)]

    # line for featureCounts
    jobstr += ['python {}/featureCounts_collector.py -p {} -g {} -o {} '.format(abspath(join_path(iLoop_RNAseq_pipeline.__path__[0], 'scripts')),
                                                                                abspath(project_path),
                                                                                abspath(join_path(join_path(project_path, 'inputs', 'groups.json'))),
                                                                                output)]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def edgeR_job(project_path, groups, output, collectjobID, defaults, ppn='1', walltime="12:00:00"):
    '''Prepare inputs for edgeR and generate a job script'''

    # Generate conditions input file
    condf = open(abspath(join_path(project_path, 'inputs', 'conditions_Rready.csv')), 'w')
    condf.writelines('Sample,Strain,Treatment\n')
    condf.writelines(
        [','.join([sample, 'strain', group + '\n']) for group, samples in groups.items() for sample in samples.keys()])
    condf.close()

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'edgeR')
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', 'afterok:{}'.format(collectjobID))
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['Rscript {}/edge_Rscript.r -p {}, -c {} -s {} -o {}'.format(abspath(join_path(iLoop_RNAseq_pipeline.__path__, 'scripts')),
                                                                           project_path,
                                                                           abspath(join_path(project_path,
                                                                                             'results',
                                                                                             'featureCounts_collected.csv')),
                                                                           abspath(join_path(project_path,
                                                                                             'inputs',
                                                                                             'conditions_Rready.csv')),
                                                                           output)]




    return '\n\n'.join(jobstr).replace('PPN', ppn)


def merge_job(project_path, mapjobIDs, ref, defaults, ppn='1', walltime='01:00:00'):
    logger.info('Using cuffmerge options: {}'.format(defaults['cuffmerge_options']))
    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'cuffmerge')
                         .replace('WALLTIME', walltime)
                         .replace('PROJECT', defaults['project'])
                         .replace('DEPEND', (
        'afterok:{}'.format(':'.join([mapjob for mapjob in mapjobIDs])) if mapjobIDs != [''] else ''))
                         .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                         .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']
    jobstr += ['cuffmerge {} {} -p PPN {} -o {} {}'.format(defaults['cuffmerge_options'],
                                                           (('-g ' + ref['gff_genome']) if ref.get(
                                                               'gff_genome') else ''),
                                                           (('-s ' + ref['fasta_genome']) if ref.get(
                                                               'fasta_genome') else ''),
                                                           (abspath(join_path(project_path, 'cmerge',
                                                                                         'merged_asm'))),
                                                           (abspath(join_path(project_path, 'cmerge',
                                                                                         'assemblies.txt'))))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def quant_jobs(project_path, sample, mergejob, ref, defaults, ppn='8', walltime='12:00:00'):
    logger.info('Using cuffquant options: {}'.format(defaults['cuffquant_options']))
    jobstr = []
    jobstr += [job_header.replace('JOBNAME', '_'.join([sample] + ['cuffquant']))
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', 'afterok:{}'.format(mergejob))
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    # make this job depend on successful completion of previous jobs: merge_job
    jobstr += ['#PBS -W depend=afterok:{}'.format(mergejob)]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']

    if ref.get('bowtie_indexes'):
        jobstr += ['export BOWTIE_INDEXES={}'.format(ref['bowtie_indexes'])]

    jobstr += ['cuffquant {} -p PPN {} -o {} {} {} {} '.format(defaults['cuffquant_options'],
                                                               ('-M ' + ref['gff_mask']) if ref.get('gff_mask') else '',
                                                               (abspath(join_path(project_path, sample))),
                                                               ('-b ' + ref['fasta_genome']) if ref.get(
                                                                   'fasta_genome') else '',
                                                               (abspath(
                                                                   join_path(project_path, 'cmerge', 'merged_asm',
                                                                                'merged.gtf'))),
                                                               (abspath(join_path(project_path, sample,
                                                                                             'accepted_hits.sorted.bam'))))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def diff_job(project_path, groups, quantjobsIDs, ppn='8', walltime='24:00:00', ref=None, defaults=None):
    logger.info('Using cuffdiff options: {}'.format(defaults['cuffdiff_options']))
    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'cuffdiff')
                   .replace('WALLTIME', walltime)
                   .replace('PROJECT', defaults['project'])
                   .replace('DEPEND', (
        'afterok:{}'.format(':'.join([qjob for qjob in quantjobsIDs])) if quantjobsIDs != [''] else ''))
                   .replace('JOB_OUTPUTS', abspath(join_path(project_path, 'job_outputs')))
                   .replace('EMAILADDRESS', defaults['email'])]

    # make this job depend on successful completion of previous jobs: mapandlink_jobs
    jobstr += ['#PBS -W depend=afterok:{}'.format(':'.join([quantjobsID for quantjobsID in quantjobsIDs]))]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']

    jobstr += ['cuffdiff {} -o {} {} -p PPN {} -L {} -u {} -o {} {}'.format(defaults['cuffdiff_options'],
                                                                            abspath(
                                                                                join_path(project_path, 'cdiff',
                                                                                             'diff_out')),
                                                                            (('-b ' + ref['fasta_genome']) if ref['fasta_genome'] != '' else ''),
                                                                            (('-M ' + ref['gff_mask']) if ref['gff_mask'] != '' else ''),
                                                                            (','.join(
                                                                                [group_name for group_name, group
                                                                                 in groups.items()])),
                                                                            (abspath(
                                                                                join_path(project_path, 'cmerge',
                                                                                             'merged_asm',
                                                                                             'merged.gtf'))),
                                                                            (abspath(
                                                                                join_path(project_path, 'cdiff',
                                                                                             'diff_out'))),
                                                                            (' '.join(
                                                                                [','.join([abspath(join_path(
                                                                                    project_path, sample,
                                                                                    'abundances.cxb')) for
                                                                                           sample, reads in
                                                                                           group.items()])
                                                                                 for group_name, group in
                                                                                 groups.items()])))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def job_submitter(js, path, name):
    jfn = join_path(path, name)
    jf = open(jfn, 'w')
    jf.write(js)
    jf.close()
    p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    os.system('sleep 0.5')
    return out.strip().decode(sys.getdefaultencoding())

def job_organizer(project_path, groups, ref, defaults, map_to_mask, ppn='8', readtype='raw', jobs=None):
    try:
        job_files_path = abspath(join_path(project_path, 'job_files'))
        os.mkdir(job_files_path)
    except FileExistsError:
        logger.warning('Folder for job files exists. Previously generated job files will be overwritten.')

    try:
        os.mkdir(join_path(project_path, 'job_outputs'))
    except FileExistsError:
        logger.warning('Folder for job outputs exists.')

    try:
        results_path = abspath(join_path(project_path, 'results'))
        os.mkdir(join_path(results_path))
    except FileExistsError:
        logger.warning('Folder for results exists.')

    # do quality checks
    if ('fastqc' in jobs) or (jobs == []):
        try:
            output_path = abspath(join_path(project_path, 'reads', 'QC_output', readtype))
            if not check_fastqc(groups=groups, output_path=output_path):
                js = fastqc_job(project_path=project_path, groups=groups, output_path=output_path, defaults=defaults)
                qcjobID = job_submitter(js, job_files_path, 'job_fastqc.sh')
            else:
                logger.info('Existing fastqc files found. Skipping quality check job.')
                qcjobID = True
        except Exception as ex:
            logger.error(
                'Problem with FastQC. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False

    # TODO: Come up with a qc threshold to continue or terminate jobs. Use qcjobID from above.

    # generate and submit map and link jobs
    if map_to_mask:
        try:
            for group_name, group in groups.items():
                for sample, reads in group.items():
                    js = mapandlink_jobs(project_path=project_path, sample=sample, reads=reads, ref=ref,
                                         defaults=defaults, ppn=ppn, jobs=jobs, map_to_mask=map_to_mask)
                    job_submitter(js, job_files_path, 'job_{}_mapandlink_to_mask.sh'.format(sample))
        except Exception as ex:
            logger.error(
                'Problem with map and link. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False
    # TODO: a threshold for mapping to mask/genome can be used to decide whether to go forward or stop pipeline. At least generate a warning.

    mljobs = ['hisat2', 'stringtie', 'cufflinks', 'htseq-count', 'featureCounts']
    if (any(job for job in jobs if job in mljobs)) or (jobs == []):
        try:
            mapjobIDs = []
            for group_name, group in groups.items():
                for sample, reads in group.items():
                    js = mapandlink_jobs(project_path=project_path, sample=sample, reads=reads, ref=ref,
                                         defaults=defaults, ppn=ppn, jobs=jobs)
                    mapjobIDs.append(job_submitter(js, job_files_path, 'job_{}_mapandlink.sh'.format(sample)))
        except Exception as ex:
            logger.error(
                'Problem with map and link. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False

        # TODO: this check should wait for the end of jobs. Ideally, should be included in their respective job files. Ignoring for now.
        # if not check_mapandlink(groups, project_path):
        #     return False

    else:
        mapjobIDs = ['']

    # collect htseq counts - Ignore this since featureCounts is prefered now
    # if any(job for job in jobs if job in ['htseq-count', 'edgeR', 'DESeq', 'htseq-count-collect']) or (jobs == []):
    #     try:
    #         js = collect_counts_job(project_path=project_path, output=abspath(join_path(results_path, 'htseq_counts_collected')), mapjobIDs=mapjobIDs, defaults=defaults)
    #         collectjobID = job_submitter(js=js, path=job_files_path, name='job_htseq_count_collector.sh')
    #     except Exception as ex:
    #         logger.error(
    #             'Problem with HTseq-count collector. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
    #                 type(ex).__name__, ex.args))
    #         return False

    if any(job for job in jobs if job in ['featureCounts', 'edgeR', 'DESeq', 'featureCounts-collect']) or (jobs == []):
        try:
            js = collect_counts_job(project_path=project_path, output=abspath(join_path(results_path, 'featureCounts_collected')), mapjobIDs=mapjobIDs, defaults=defaults)
            collectjobID = job_submitter(js=js, path=job_files_path, name='job_featureCounts_collector.sh')
        except Exception as ex:
            logger.error(
                'Problem with featureCounts collector. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False

    if 'edgeR' in jobs:
        try:
            js = edgeR_job(project_path=project_path, groups=groups, output=results_path, collectjobID=collectjobID, defaults=defaults)
            edgeRjobID= job_submitter(js=js, path=job_files_path, name='job_edgeR.sh')
        except Exception as ex:
            logger.error(
                'Problem with edgeR. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False

    # generate and submit merge job
    if ('cuffmerge' in jobs) or (jobs == []):
        try:
            os.mkdir(join_path(project_path, 'cmerge'))
        except FileExistsError:
            logger.warning('Folder for cuffmerge exists. Previously generated files will be overwritten.')

        try:
            af = open(join_path(project_path, 'cmerge', 'assemblies.txt'), 'w')
            af.write('\n'.join([abspath(join_path(project_path, replicate, 'transcripts.gtf')) for group in groups.values()
                                for replicate in group.keys()]))
            af.close()
            logger.info('"Assemblies.txt" is generated under "cmerge".')
        except Exception as ex:
            logger.error(
                '"Assemblies.txt for cuffmerge could not be generated.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False

        try:
            js = merge_job(project_path=project_path, mapjobIDs=mapjobIDs, ref=ref, defaults=defaults)
            mergejob = job_submitter(js, job_files_path, 'job_cuffmerge.sh')
        except Exception as ex:
            logger.error(
                'Problem with Cuffmerge. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False
    else:
        mergejob = ''

    # generate and submit cuffquant jobs
    if ('cuffquant' in jobs) or (jobs == []):
        try:
            for group_name, group in groups.items():
                quantjobsIDs = []
                for sample, reads in group.items():
                    js = quant_jobs(project_path=project_path, sample=sample, mergejob=mergejob, ref=ref,
                                    defaults=defaults, ppn=ppn)
                    quantjobsIDs.append(job_submitter(js, job_files_path, 'job_{}_cuffquant.sh'.format(sample)))
        except Exception as ex:
            logger.error(
                'Problem with Cuffquant jobs. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False
    else:
        quantjobsIDs = ['']

    # generate and submit cuffdiff job
    if ('cuffdiff' in jobs) or (jobs == []):
        try:
            os.mkdir(join_path(project_path, 'cdiff'))
        except FileExistsError:
            logger.warning('Folder for cuffdiff exists. Previously generated files will be overwritten.')

        try:
            js = diff_job(project_path=project_path, groups=groups, quantjobsIDs=quantjobsIDs, ppn=ppn,
                          walltime='24:00:00',
                          ref=ref, defaults=defaults)
            diffjob = job_submitter(js, job_files_path, 'job_cuffdiff.sh')
        except Exception as ex:
            logger.error(
                'Problem with Cuffdiff. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(
                    type(ex).__name__, ex.args))
            return False


    logger.info('All jobs are completed.')
    return True
