#!/usr/bin/env python

import argparse
import logging
import os
import sys
import sqlalchemy
from cdis_pipe_utils import pipe_util

import tools.mutect_tool as mutect_tool


def is_nat(x):
    '''
    Checks that a value is a natural number.
    '''
    if int(x) > 0:
        return int(x)
    raise argparse.ArgumentTypeError('%s must be positive, non-zero' % x)

def main():
    parser = argparse.ArgumentParser('GATK MuTect2 Panel Of Normal creation')

    # Logging flags.
    parser.add_argument('-d', '--debug',
        action = 'store_const',
        const = logging.DEBUG,
        dest = 'level',
        help = 'Enable debug logging.',
    )
    parser.set_defaults(level = logging.INFO)

    # Required flags.

    parser.add_argument('-r', '--reference_fasta_path',
                        required = False,
                        help = 'Reference fasta path.',
    )

    parser.add_argument('-rf', '--reference_fasta_fai',
                        required = False,
                        help = 'Reference fasta fai path.',
    )

    parser.add_argument('-snp', '--known_snp_vcf_path',
                        required = False,
                        help='Reference SNP path.',
    )

    parser.add_argument('-cos', '--cosmic_path',
                        required = False,
                        help='Reference COSMIC path.',
    )

    parser.add_argument('-n', '--normal_bam_path',
                        required = False,
                        help = 'normal bam path.',
    )

    parser.add_argument('-t', '--tumor_bam_path',
                        required = False,
                        help = 'tumor bam path',
    )

    parser.add_argument('-v', '--vcf_path',
                        required = False,
                        help = 'Individual VCF path',
    )

    parser.add_argument('-j', '--thread_count',
                        required = False,
                        type = is_nat,
                        help = 'Maximum number of threads for execution.',
    )

    parser.add_argument('-bs', '--Parallel_Block_Size',
                        type = is_nat,
                        default = 50000000,
                        required = False,
                        help = 'Parallel Block Size',
    )

    parser.add_argument('-u', '--uuid',
                        required = True,
                        help = 'analysis_id string',
    )

    parser.add_argument('--normal_panel',
                        required = True,
                        help = 'panel of normals vcf'
    )

    parser.add_argument('--contEst',
                        required = True,
                        help = 'Contamination estimation value from ContEst'
    )

    parser.add_argument('--case_id',
                        required = True,
                        help = 'case id for the tumor-normal pair'
    )

    args = parser.parse_args()
    uuid = args.uuid
    thread_count = str(args.thread_count)
    Parallel_Block_Size = str(args.Parallel_Block_Size)

    logger = pipe_util.setup_logging('gatk_mutect2', args, uuid)
    engine = pipe_util.setup_db(uuid)

    hostname = os.uname()[1]
    logger.info('hostname=%s' % hostname)


    normal_bam_path = pipe_util.get_param(args, 'normal_bam_path')
    tumor_bam_path = pipe_util.get_param(args, 'tumor_bam_path')
    known_snp_vcf_path = pipe_util.get_param(args, 'known_snp_vcf_path')
    cosmic_path = pipe_util.get_param(args, 'cosmic_path')
    reference_fasta_path = pipe_util.get_param(args, 'reference_fasta_path')
    thread_count = pipe_util.get_param(args, 'thread_count')
    fai_path = pipe_util.get_param(args, 'reference_fasta_fai')
    blocksize = pipe_util.get_param(args, 'Parallel_Block_Size')
    normal_panel = pipe_util.get_param(args, 'normal_panel')
    contEst = pipe_util.get_param(args, 'contEst')

    mutect_tool.run_mutect(uuid,
                       normal_bam_path,
                       tumor_bam_path,
                       normal_panel,
                       contEst,
                       thread_count,
                       reference_fasta_path,
                       cosmic_path,
                       fai_path,
                       blocksize,
                       known_snp_vcf_path,
                       engine,
                       logger)


if __name__ == '__main__':

    main()
