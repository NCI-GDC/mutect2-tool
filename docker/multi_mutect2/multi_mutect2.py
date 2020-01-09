'''
Multithreading MuTect2

@author: Shenglai Li
'''

import os
import sys
import time
import logging
import argparse
import subprocess
import string
from functools import partial
from multiprocessing.dummy import Lock, Pool

def setup_logger():
    '''
    Sets up the logger.
    '''
    logger = logging.getLogger("multi_mutect2")
    logger_format = '[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s'
    logger.setLevel(level=logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(logger_format, datefmt='%Y%m%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def do_pool_commands(cmd, logger, shell_var=True, lock=Lock()):
    '''run pool commands'''
    try:
        output = subprocess.Popen(
            cmd,
            shell=shell_var,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        output_stdout, output_stderr = output.communicate()
        with lock:
            logger.info('MuTect2 Args: %s', cmd)
            logger.info(output_stdout)
            logger.info(output_stderr)
    except BaseException:
        logger.error("command failed %s", cmd)
    return output.wait()

def multi_commands(cmds, thread_count, logger, shell_var=True):
    '''run commands on number of threads'''
    pool = Pool(int(thread_count))
    output = pool.map(partial(do_pool_commands, logger=logger, shell_var=shell_var), cmds)
    return output

def get_region(intervals):
    '''get region from intervals'''
    interval_list = []
    with open(intervals, 'r') as fh:
        line = fh.readlines()
        for bed in line:
            blocks = bed.rstrip().rsplit('\t')
            intv = '{}:{}-{}'.format(blocks[0], int(blocks[1])+1, blocks[2])
            interval_list.append(intv)
    return interval_list

def cmd_template(dct):
    '''cmd template'''
    lst = [
        'java', '-Djava.io.tmpdir=/tmp/job_tmp_${BLOCK_NUM}',
        '-d64', '-jar', '-Xmx${JAVA_HEAP}', '-XX:+UseSerialGC',
        '${GATK_PATH}', '-T', 'MuTect2', '-nct', '1', '-nt', '1',
        '-R', '${REF}',
        '-L', '${REGION}',
        '-I:tumor', '${TUMOR_BAM}',
        '-I:normal', '${NORMAL_BAM}',
        '--normal_panel', '${PON}',
        '--cosmic', '${COSMIC}',
        '--dbsnp', '${DBSNP}',
        '--contamination_fraction_to_filter', '${CONTAMINATION}',
        '-o', '${BLOCK_NUM}.mt2.vcf',
        '--output_mode', 'EMIT_VARIANTS_ONLY',
        '--disable_auto_index_creation_and_locking_when_reading_rods'
    ]
    if dct['dontUseSoftClippedBases']:
        lst = lst + ['--dontUseSoftClippedBases']
    template = string.Template(' '.join(lst))
    for i, interval in enumerate(get_region(dct['interval_bed_path'])):
        cmd = template.substitute(
            dict(
                BLOCK_NUM=i,
                JAVA_HEAP=dct['java_heap'],
                GATK_PATH=dct['gatk_path'],
                REF=dct['reference_path'],
                REGION=interval,
                TUMOR_BAM=dct['tumor_bam'],
                NORMAL_BAM=dct['normal_bam'],
                PON=dct['pon'],
                COSMIC=dct['cosmic'],
                DBSNP=dct['dbsnp'],
                CONTAMINATION=dct['contest']
            )
        )
        yield cmd, '{}.mt2.vcf'.format(i)

def get_args():
    '''
    Loads the parser.
    '''
    # Main parser
    parser = argparse.ArgumentParser('Internal multithreading MuTect2 calling.')
    # Required flags.
    parser.add_argument(
        '-j',
        '--java_heap',
        required=True,
        help='Java heap memory.'
    )
    parser.add_argument(
        '-f',
        '--reference_path',
        required=True,
        help='Reference path.'
    )
    parser.add_argument(
        '-r',
        '--interval_bed_path',
        required=True,
        help='Interval bed file.'
    )
    parser.add_argument(
        '-t',
        '--tumor_bam',
        required=True,
        help='Tumor bam file.'
    )
    parser.add_argument(
        '-n',
        '--normal_bam',
        required=True,
        help='Normal bam file.'
    )
    parser.add_argument(
        '-c',
        '--thread_count',
        type=int,
        required=True, help='Number of thread.'
    )
    parser.add_argument(
        '-p',
        '--pon',
        required=True,
        help='Panel of normals reference path.'
    )
    parser.add_argument(
        '-s',
        '--cosmic',
        required=True,
        help='Cosmic reference path.'
    )
    parser.add_argument(
        '-d',
        '--dbsnp',
        required=True,
        help='dbSNP reference path.'
    )
    parser.add_argument(
        '-e',
        '--contest',
        required=True,
        help='Contamination estimation value from ContEst.'
    )
    parser.add_argument(
        '-m',
        '--dontUseSoftClippedBases',
        action="store_true"
        , help='If specified, it will not analyze soft clipped bases in the reads.'
    )
    return parser.parse_args()

def main(args, logger):
    '''main'''
    logger.info("Running GATK3.6 MuTect2")
    dct = vars(args)
    dct['gatk_path'] = '/opt/GenomeAnalysisTK.jar'
    mutect2_cmds = list(cmd_template(dct))
    outputs = multi_commands([x[0] for x in mutect2_cmds], dct['thread_count'], logger)
    if any(x != 0 for x in outputs):
        logger.error('Failed multi_mutect2')
    else:
        logger.info('Completed multi_mutect2')
        first = True
        with open('multi_mutect2_merged.vcf', 'w') as oh:
            for _, out in mutect2_cmds:
                with open(out) as fh:
                    for line in fh:
                        if first or not line.startswith('#'):
                            oh.write(line)
                first = False

if __name__ == '__main__':
    # CLI Entrypoint.
    start = time.time()
    logger_ = setup_logger()
    logger_.info('-'*80)
    logger_.info('multi_mutect2.py')
    logger_.info('Program Args: %s', ' '.join(sys.argv))
    logger_.info('-'*80)

    args_ = get_args()

    # Process
    logger_.info(
        'Processing tumor and normal bam files %s, %s',
        os.path.basename(args_.tumor_bam),
        os.path.basename(args_.normal_bam)
    )
    main(args_, logger_)

    # Done
    logger_.info('Finished, took %s seconds', round(time.time() - start, 2))
