from __future__ import division, print_function
import os
import subprocess
import logging

logger = logging.getLogger(__name__)

job_header = \
'''#!/bin/sh
### Account information
#PBS -W group_list=PROJECT -A PROJECT
# -- Name of the job ---
#PBS -N JOBNAME
# -- estimated wall clock time (execution time): hh:mm:ss --
#PBS -l walltime=WALTIME
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
        if (not os.path.isfile(os.path.join(output_path, qcf + '_fastqc.zip'))) or \
           (not os.path.isfile(os.path.join(output_path, qcf + '_fastqc.html'))):
            return False
        elif (os.path.getsize(os.path.join(output_path, qcf + '_fastqc.zip')) == 0) or \
             (os.path.getsize(os.path.join(output_path, qcf + '_fastqc.html')) == 0):
            return False
    return True

def fastqc_job(project_path, groups, output_path, defaults, ppn='8', walltime ='02:00:00', ):
    """Runs fastqc for all the reads.

    reads: dictionary from set_project function
    output_path: path for the output of fastqc
                - ideally, project_path/reads/QC_output/XYZ where XYZ denotes raw or trimmed
     """
    # runs a quality check on fastq files

    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'fastqc') \
                   .replace('WALTIME', walltime) \
                   # .replace('PPN', ppn)\
                   .replace('PROJECT', defaults['project']) \
                   .replace('JOB_OUTPUTS',  os.path.join(project_path, 'job_outputs')) \
                   .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs FastQC/0.11.2''']

    jobstr += ['fastqc -t PPN -o {} {}'.format(output_path,
                                               ' '.join([readf for group in groups.values() for readfs in group.values() for readf in readfs]))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def mapandlink_jobs(project_path, sample, reads, defaults, ref, ppn='8', walltime ='12:00:00', jobs=None):

    if jobs is None:
        jobs = ['hisat2', 'stringtie', 'cufflinks', 'htseq-count']

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', '_'.join([sample]+jobs))\
                       .replace('WALTIME', walltime)\
                       # .replace('PPN', ppn)\
                       .replace('PROJECT', defaults['project']) \
                       .replace('JOB_OUTPUTS',  os.path.join(project_path, 'job_outputs')) \
                       .replace('EMAILADDRESS', defaults['email'])]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools samtools/1.2 bowtie2/2.2.5 cufflinks/2.2.1
# hisat2 and stringtie modules are not kept up-to-date, so use the local ones:
PATH=$PATH:/home/projects/cu_10010/programs/hisat2-2.0.4:/home/projects/cu_10010/programs/stringtie-1.2.3.Linux_x86_64
export PATH''']

    jobstr += ['export HISAT2_INDEXES={}'.format(ref['hisat2_indexes'])]
    jobstr += ['mkdir {}/{}'.format(project_path, sample)]

    R1reads = [read for read in reads if read[-15:-13] == 'R1']
    R1reads.sort()
    R2reads = [read for read in reads if read[-15:-13] == 'R2']
    R2reads.sort()

    if (ref.get('hisat2_indexes')) and ('hisat2' in jobs):
        jobstr += ['''echo "hisat2"
hisat2 {} -p PPN -x {} {} {} 2>{} | \
samtools view -@ PPN -hbu - | \
samtools sort -@ PPN - {}'''.format(defaults['hisat2_options'],
                                    ref['hisat2_indexes'],
                                    ('-1'+','.join(R1reads)),
                                    ('-2'+','.join(R2reads)),
                                    (os.path.join(project_path, sample, 'align_summary.txt')),
                                    (os.path.join(project_path, sample, 'accepted_hits.sorted')))]

    if 'stringtie' in jobs:
        logger.warning('Beware: Stringtie does not allow masking for now.')
        jobstr += ['echo "stringtie"\nstringtie {} -p PPN {} -o {} -A {} {}'.format(defaults['stringtie_options'],
                                                                  (('-G ' + ref['gff_genome']) if ref.get('gff_genome') else ''),
                                                                  os.path.join(project_path, sample, 'transcripts.gtf'),
                                                                  os.path.join(project_path, sample, 'gene_abund.tab'),
                                                                  os.path.join(project_path, sample, 'accepted_hits.sorted.bam'))]

    if 'cufflinks' in jobs:
        jobstr += ['echo "cufflinks"\ncufflinks {} -p PPN {} {} -o {} {}'.format(defaults['cufflinks_options'],
                                                                ('-G '+ref['gff_genome']) if ref.get('gff_genome') else '',
                                                                ('-M ' + ref['gff_mask']) if ref.get('gff_mask') else '',
                                                                (os.path.join(project_path, sample)),
                                                                (os.path.join(project_path, sample, 'accepted_hits.sorted.bam')))]

    if 'htseq-count' in jobs:
        jobstr += ['echo "htseq"\nhtseq-count {} -f bam {} {} -o {} > {}'.format(defaults['htseq_options'],
                                                                   (os.path.join(project_path, sample, 'accepted_hits.sorted.bam')),
                                                                   ((ref['gff_genome']) if ref.get('gff_genome') else ''),
                                                                   (os.path.join(project_path, sample, 'htseq_counts.sam')),
                                                                   (os.path.join(project_path, sample, 'htseq_counts.out')))]

    return '\n\n'.join(jobstr).replace('PPN', str(ppn))


def merge_job(project_path, mapjobIDs, ppn='1', walltime='01:00:00', ref=None, defaults=None):

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'cuffmerge')\
                         .replace('WALTIME', walltime)\
                         # .replace('PPN', ppn)\
                         .replace('PROJECT', defaults['project']) \
                         .replace('JOB_OUTPUTS',  os.path.join(project_path, 'job_outputs')) \
                         .replace('EMAILADDRESS', defaults['email'])]

    # make this job depend on successful completion of previous jobs: mapandlink_jobs
    jobstr += ['#PBS –W depend=afterok:{}'.format(':'.join([mapjob for mapjob in mapjobIDs]))]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']

    jobstr += ['cuffmerge {} {} -p PPN -o {} assemblies.txt'.format(defaults['cufflinks_options'],
                                                                    (('-g '+ref['gff_genome']) if ref.get('gff_genome') else ''),
                                                                    (('-s ' + ref['fasta_genome']) if ref.get('fasta_genome') else ''),
                                                                    (os.path.join(project_path, 'cmerge', 'merged_asm')))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def quant_jobs(project_path, sample, mergejob, ppn='8', walltime ='12:00:00', ref=None, defaults=None):

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', '_'.join([sample]+'cuffquant'))\
                       .replace('WALTIME', walltime)\
                       # .replace('PPN', ppn)\
                       .replace('PROJECT', defaults['project']) \
                       .replace('JOB_OUTPUTS',  os.path.join(project_path, 'job_outputs')) \
                       .replace('EMAILADDRESS', defaults['email'])]

    # make this job depend on successful completion of previous jobs: merge_job
    jobstr += ['#PBS –W depend=afterok:{}'.format(mergejob)]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']

    if ref.get('bowtie_indexes'):
        jobstr += ['export BOWTIE_INDEXES={}}'.format(ref['bowtie_indexes'])]

    jobstr += ['cuffquant {} -p PPN {} -o {} {} {} '.format(defaults['cuffquant_options'],
                                                            ('-M ' + ref['gff_mask']) if ref.get('gff_mask') else '',
                                                            (os.path.join(project_path, sample, 'accepted_hits.sorted.bam')),
                                                            ('-b ' + ref['fasta_genome']) if ref.get('fasta_genome') else '',
                                                            (os.path.join(project_path, 'cmerge', 'merged_asm', 'merged.gtf')),
                                                            (os.path.join(project_path, sample, 'accepted_hits.sorted.bam')))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def diff_job(project_path, groups, quantjobsIDs, ppn='8', walltime='24:00:00', ref=None, defaults=None):

    jobstr = []
    jobstr += [job_header.replace('JOBNAME', 'cuffdiff')\
                         .replace('WALTIME', walltime)\
                         # .replace('PPN', ppn)\
                         .replace('PROJECT', defaults['project']) \
                         .replace('JOB_OUTPUTS',  os.path.join(project_path, 'job_outputs')) \
                         .replace('EMAILADDRESS', defaults['email'])]

    # make this job depend on successful completion of previous jobs: mapandlink_jobs
    jobstr += ['#PBS –W depend=afterok:{}'.format(':'.join([quantjobsID for quantjobsID in quantjobsIDs]))]

    jobstr += ['''# Load modules needed by myapplication.x
module load ngs tools cufflinks/2.2.1 tophat/2.1.1 bowtie2/2.2.5''']

    jobstr += ['cuffdiff {} -o diff_out {} -p PPN {} -L {} -u {} -o {} {}'.format(ref['cuffdiff_options'],
                                                                                  (('-b '+ref['fasta_genome']) if ref['fasta_genome'] != '' else ''),
                                                                                  (('-M '+ref['gff_mask']) if ref['gff_mask'] != '' else ''),
                                                                                  (','.join([group_name for group_name, group in groups.items()])),
                                                                                  (os.path.join(project_path, 'cmerge', 'merged_asm', 'merged.gtf')),
                                                                                  (os.path.join(project_path, 'cdiff', 'diff_out')),
                                                                                  (' '.join([','.join([os.path.join(project_path, sample, 'abundances.cxb') for sample, reads in group.items()]) for group_name, group in groups.items()])))]

    return '\n\n'.join(jobstr).replace('PPN', ppn)


def job_submitter(project_path, groups, ref, defaults, ppn='8', readtype='raw'):

    try:
        os.mkdir(os.path.join(project_path, 'job_files'))
    except FileExistsError:
        logger.warning('Folder for job files exists. Previously generated job files will be overwritten.')

    try:
        os.mkdir(os.path.join(project_path, 'job_outputs'))
    except FileExistsError:
        logger.warning('Folder for job outputs exists.')

    # do quality checks
    try:
        output_path = os.path.join(project_path, 'reads', 'QC_output', readtype)
        if not check_fastqc(groups=groups, output_path=output_path):
            js = fastqc_job(project_path=project_path, groups=groups, output_path=output_path, defaults=defaults)
            jfn = os.path.join(project_path, 'job_files', 'job_fastqc.sh')
            jf = open(jfn, 'w')
            jf.write(js)
            jf.close()
            p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
            out, err = p.communicate()
            qcjobID = out.split(b'.')[0]
            os.system('sleep 1')
        else:
            logger.info('Existing fastqc files found. Skipping quality check job.')
            qcjobID = True
    except Exception as ex:
        logger.error('Problem with FastQC. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(type(ex).__name__, ex.args))
        return False

    # TODO: Come up with a qc threshold to continue or terminate jobs. Use qcjobID from above.

    # generate and submit map and link jobs
    try:
        for group_name, group in groups.items():
            mapjobIDs = []
            samples = []
            for sample, reads in group.items():
                samples.append(sample)
                js = mapandlink_jobs(project_path=project_path, sample=sample, reads=reads, ref=ref,
                                     defaults=defaults, ppn=ppn)
                jfn = os.path.join(project_path, 'job_files', 'job_{}_mapandlink.sh'.format(sample))
                jf = open(jfn, 'w')
                jf.write(js)
                jf.close()
                p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p.wait()
                out, err = p.communicate()
                mapjobIDs.append(out.split(b'.')[0])
                os.system('sleep 1')
    except Exception as ex:
        logger.error('Problem with map and link. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(type(ex).__name__, ex.args))
        return False

    # generate and submit merge job
    try:
        os.mkdir(os.path.join(project_path, 'cmerge'))
    except FileExistsError:
        logger.warning('Folder for cuffmerge exists. Previously generated files will be overwritten.')

    try:
        af = open(os.path.join([project_path, 'cmerge', 'assemblies.txt']), 'w')
        af.write('\n'.join([os.path.join(project_path, sample, 'transcripts.gtf') for sample in samples]))
        af.close()
    except:
        return False

    try:
        js = merge_job(project_path=project_path, mapjobs=mapjobIDs, ref=ref, defaults=defaults)
        jfn = os.path.join(project_path, 'job_files', 'job_cuffmerge.sh')
        jf = open(jfn, 'w')
        jf.write(js)
        jf.close()
        p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
        out, err = p.communicate()
        mergejob = out.split(b'.')[0]
        os.system('sleep 1')
    except Exception as ex:
        logger.error('Problem with Cuffmerge. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(type(ex).__name__, ex.args))
        return False

    # generate and submit cuffquant jobs
    try:
        for group_name, group in groups.items():
            quantjobsIDs = []
            for sample, reads in group.items():
                js = quant_jobs(project_path=project_path, sample=sample, mergejob=mergejob, ref=ref, defaults=defaults, ppn=ppn)
                jfn = os.path.join(project_path, 'job_files', 'job_{}_cuffquant.sh'.format(sample))
                jf = open(jfn, 'w')
                jf.write(js)
                jf.close()
                p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p.wait()
                out, err = p.communicate()
                quantjobsIDs.append(out.split(b'.')[0])
                os.system('sleep 1')
    except Exception as ex:
        logger.error('Problem with Cuffquant jobs. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(type(ex).__name__, ex.args))
        return False

    # generate and submit cuffdiff job
    try:
        os.mkdir(os.path.join(project_path, 'cdiff'))
    except FileExistsError:
        logger.warning('Folder for cuffdiff exists. Previously generated files will be overwritten.')

    try:
        js = diff_job(project_path=project_path, groups=groups, quantjobsIDs=quantjobsIDs, ppn=ppn, walltime='24:00:00',
                      ref=ref, defaults=defaults)
        jfn = os.path.join(project_path, 'job_files', 'job_cuffdiff.sh')
        jf = open(jfn, 'w')
        jf.write(js)
        jf.close()
        p = subprocess.Popen(['qsub', jfn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
        out, err = p.communicate()
        diffjob = out.split(b'.')[0]
        os.system('sleep 1')
    except Exception as ex:
        logger.error('Problem with Cuffdiff. RNAseq analysis is stopped.\nAn exception of type {} occured. Arguments:\n{}'.format(type(ex).__name__, ex.args))
        return False

    return True



#
#
# def run_fastx_trimmer(name, input_path, output_path, step):
#     # trims 10 first bases of each read, returns files with "trimmed" prefix in chosen folder
#
#     data_files = glob.glob(os.path.join(input_path, '*.gz'))
#     task_count = len(data_files)
#     assert task_count > 0, "Could not find any .gz files in folder %s" % input_path
#
#     output_path = os.path.join(output_path, '2.FASTX')
#     create_path_if_not_exists(output_path)
#
#     mem_req = "10G"
#     time_req = "10:00:00"
#
#     command = "gzip -cd $input | fastx_trimmer -f10 -Q33 | gzip -c > $OUT/trimmed_${input##*/}"
#
#     write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)
#
#
# def run_star(name, input_path, output_path, step):
#     paired_q = raw_input("Do you have paired reads in separate files? [y/n] ")
#
#     if paired_q.lower() == 'n':
#         print('Ok, one file per task.')
#
#         data_files = glob.glob(os.path.join(input_path, '*.fastq.gz'))
#
#     elif paired_q.lower() == 'y':
#         print ('Ok, looking for read mates. Mates must be labeled *_1.fastq.gz and *_2.fastq.gz for pairing.')
#
#         # build forward and reverse file lists, match pairs
#         forward = []
#         reverse = []
#         for read_file in glob.glob(os.path.join(input_path, '*_2.fastq.gz')):
#             reverse.append(read_file)
#         for read_file in glob.glob(os.path.join(input_path, '*_1.fastq.gz')):
#             forward.append(read_file)
#
#         pairs = []
#         for forward_file in forward:
#             for reverse_file in reverse:
#                 for_path, ext_for = forward_file.rsplit("_", 1)
#                 rev_path, ext_rev = reverse_file.rsplit("_", 1)
#                 if for_path == rev_path:
#                     pairs.append((forward_file, reverse_file))
#
#         # double check correct pairing in  easy to read table
#         easy_read = []
#         for i in range(len(pairs)):
#             path1, name1 = pairs[i][0].rsplit("/", 1)
#             path2, name2 = pairs[i][1].rsplit("/", 1)
#             easy_read.append((i+1, name1, name2))
#
#         table = tabulate(easy_read, headers=["Pair #", "Forward", "Reverse"], tablefmt="grid")
#         print (table)
#         mates_q = raw_input("Are these pairings correct? [y/n] ")
#         if mates_q.lower() != "y":
#             print("Run aborted.")
#             return
#         else:
#             filename = os.path.join(output_path, 'table_pairs.txt')
#             f = open(filename, 'w')
#             f.write(table)
#             f.close()
#
#         data_paired = []
#         for i in range(len(pairs)):
#             data_paired.append(str.join(' ', pairs[i]))
#         data_files = data_paired
#
#     else:
#         print("Run aborted.")
#         return
#
#     task_count = len(data_files)
#     assert task_count > 0, "Could not find any fastq files in folder %s" % input_path
#
#     command = 'input1=${inputs[$SGE_TASK_ID+$SGE_TASK_ID-1]}\n' \
#     'input2=${inputs[$SGE_TASK_ID+$SGE_TASK_ID]} \n' \
#     'echo "Actual input for this task is:" $input1 $input2 \n' \
#     'STAR -runThread %somethings --genomeDir /netapp/home/dreuxj/Annotation/GRCh38_Gencode24/ --readFilesIn $input1 $input2\
#     --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonical\
#     --outFilterType BySJout --outFileNamePrefix $OUT/$SGE_TASK_ID'
#
#     output_path = os.path.join(output_path, '3.STAR')
#     create_path_if_not_exists(output_path)
#     mem_req = "35G"
#     time_req = "99:00:00"
#     write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)
#
#
# def run_samtools(name, input_path, output_path, step):
#     # samtools is installed on the cluster for all users but it's an older version
#     # - need to ref whole path for my more recent version
#
#     output_path = os.path.join(output_path, '4.SAMTOOLS')
#     create_path_if_not_exists(output_path)
#
#     time_req = "48:00:00"
#     mem_req = "10G"
#
#     if step == "samtools_view":
#
#         data_files = glob.glob(os.path.join(input_path, '*.sam'))
#         command = "/netapp/home/dreuxj/bin/samtools view -b $input -o $OUT/aligned.bam "
#
#     elif step == "samtools_sort":
#
#         data_files = glob.glob(os.path.join(input_path, '*.bam'))
#         command = "/netapp/home/dreuxj/bin/samtools sort -o $OUT/sorted.aligned.bam $input"
#
#     elif step == "samtools_idx":
#
#         data_files = glob.glob(os.path.join(input_path, '*.bam'))
#         command = "/netapp/home/dreuxj/bin/samtools index $input"
#
#     elif step == "samtools_stats":
#
#         data_files = glob.glob(os.path.join(input_path, '*.bam'))
#         command = "/netapp/home/dreuxj/bin/samtools flagstat $input > ${input}_flagstats;\
#          /netapp/home/dreuxj/bin/samtools idxstats $input > ${input}_idxstats"
#
#     else:
#         return
#
#     task_count = len(data_files)
#     assert task_count > 0, "Could not find any bam files in folder %s" % input_path
#     write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)
#
#
# def run_htseq(name, input_path, output_path, step):
#     # here think carefully about strand settings - use IGV if needed (not on all reads)
#
#     data_files = glob.glob(os.path.join(input_path, '*.sorted.aligned.bam'))
#     task_count = len(data_files)
#     assert task_count > 0, "Could not find any sorted bam files in folder %s" % input_path
#
#     command = 'python -m HTSeq.scripts.count -f bam -r pos -i gene_name -q -s no $input\
#      /netapp/home/dreuxj/Annotation/GRCh38_Gencode24/gencode.v24.primary_assembly.annotation.gtf'
#
#     output_path = os.path.join(output_path, '5.HTSEQ')
#     create_path_if_not_exists(output_path)
#
#     mem_req = "5G"
#     time_req = "99:00:00"
#     write_bash_script(name, data_files, output_path, mem_req, time_req, task_count, command, step)
#
#
# def main(argv=None):
#     """Program wrapper
#     :param argv:
#     """
#     if argv is None:
#         argv = sys.argv[1:]
#     parser = create_parser()
#     args = parser.parse_args(argv)
#
#     if args.verbose:
#         LOG.setLevel(logger.INFO)
#     if args.quiet:
#         LOG.setLevel(logger.CRITICAL)
#     if args.debug:
#         LOG.setLevel(logger.DEBUG)
#
#     name = args.name
#     step = args.step
#     input_path = os.path.join(args.base_path, args.input_path)
#     output_path = os.path.join(args.base_path, args.output_path)
#
#     create_path_if_not_exists(output_path)
#
#     if step == 'fastqc':
#         run_fastqc(name, input_path, output_path, step)
#     elif step == 'fastx':
#         run_fastx_trimmer(name, input_path, output_path, step)
#     elif step == 'star':
#         run_star(name, input_path, output_path, step)
#     elif step == 'samtools_view' or step == 'samtools_sort' or step == 'samtools_idx' or step == 'samtools_stats':
#         run_samtools(name, input_path, output_path, step)
#     elif step == 'htseq':
#         run_htseq(name, input_path, output_path, step)
#
#     else:
#         LOG.error('Did not understand step "%s". Possible values are as follows:\
#          fastqc, fastx, star, htseq,samtools_view, samtools_sort, samtools_idx, samtools_stats, Run aborted.' % step)
#         return 1
#     return 0
#
# if __name__ == '__main__':
#     import doctest
#     doctest.testmod()
#     sys.exit(main())
