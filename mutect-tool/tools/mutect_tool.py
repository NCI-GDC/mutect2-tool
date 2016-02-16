import os
import sys
import string
import subprocess
from itertools import islice
from functools import partial
from multiprocessing.dummy import Pool, Lock
from cdis_pipe_utils import df_util
from cdis_pipe_utils import pipe_util
from cdis_pipe_utils import time_util

def do_pool_commands(cmd, uuid, engine, logger, lock = Lock()):
    logger.info('running mutect2_pon chunk call: %s' % cmd)
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output_stdout = output.communicate()[1]
    with lock:
        logger.info('contents of output=%s' % output_stdout.decode().format())
        df = time_util.store_time(uuid, cmd, output_stdout, logger)
        df['cmd'] = cmd
        unique_key_dict = {'uuid': uuid, 'cmd': cmd}
        table_name = 'time_mem_MuTect2_Panel_Of_Normal_chunk_call_processes'
        df_util.save_df_to_sqlalchemy(df, unique_key_dict, table_name, engine, logger)
        logger.info('completed mutect2_pon chunk call: %s' % str(cmd))
    return output.wait()

def multi_commands(uuid, cmds, thread_count, engine, logger):
    pool = Pool(int(thread_count))
    output = pool.map(partial(do_pool_commands, uuid=uuid, engine=engine, logger=logger), cmds)
    return output

def fai_chunk(fai_path, blocksize):
  seq_map = {}
  with open(fai_path) as handle:
    head = list(islice(handle, 25))
    for line in head:
      tmp = line.split("\t")
      seq_map[tmp[0]] = int(tmp[1])
    for seq in seq_map:
        l = seq_map[seq]
        for i in range(1, l, blocksize):
            yield (seq, i, min(i+blocksize-1, l))

def mutect2_cmd_template(gatk_path, ref, fai_path, blocksize, cosmic,
                         dbsnp, output_base, normal_panel, contEst,
                         normal_bam_path, tumor_bam_path):

    template = string.Template("/usr/bin/time -v java -Djava.io.tmpdir=/tmp/job_tmp -d64 -jar ${GATK_PATH} -T MuTect2 -R ${REF} -L ${REGION} -I:tumor ${TUMOR_BAM} -I:normal ${NORMAL_BAM} --cosmic ${COSMIC} --dbsnp ${DBSNP} --artifact_detection_mode --contamination_fraction_to_filter ${CONTAMINATION} --normal_panel ${PON} -o ${OUTPUT_BASE}.${BLOCK_NUM}.mt2pon.vcf")

    for i, block in enumerate(fai_chunk(fai_path, blocksize)):

        cmd = template.substitute(
                              dict(
                                   REF = ref,
                                   REGION = '%s:%s-%s' % (block[0], block[1], block[2]),
                                   BLOCK_NUM = i),
                                   GATK_PATH = gatk_path,
                                   NORMAL_BAM = normal_bam_path,
                                   TUMOR_BAM = tumor_bam_path,
                                   CONTAMINATION = contEst,
                                   PON = normal_panel,
                                   COSMIC = cosmic,
                                   DBSNP = dbsnp,
                                   OUTPUT_BASE = output_base
        )

        yield cmd

def run_mutect(uuid, normal_bam_path, tumor_bam_path, normal_panel,
              contEst, thread_count, reference_fasta_path,cosmic_path,
              fai_path, blocksize, known_snp_vcf_path, engine, logger):

    step_dir = os.path.join(os.getcwd(), 'mutect')

    if not os.path.isdir(step_dir):
        os.makedirs(step_dir)

    normal_bam_name = os.path.basename(normal_bam_path)
    normal_base, normal_ext = os.path.splitext(normal_bam_name)

    tumor_bam_name = os.path.basename(tumor_bam_path)
    tumor_base, tumor_ext = os.path.splitext(tumor_bam_name)

    merge_dir = os.getcwd()

    #os.makedirs(merge_dir, exist_ok=True)

    #out_vcf = os.path.join(merge_dir, cram_base) + '_pon.vcf'
    #logger.info('out_vcf=%s' % out_vcf)

    #Run MuTect
    logger.info('running MuTect for: %s and %s' %(normal_bam_path, tumor_bam_path))
    home_dir = os.path.expanduser('~')
    gatk_path = os.path.join(home_dir, 'tools', 'GenomeAnalysisTK.jar')

    cmds = list(mutect2_cmd_template(
                                   gatk_path = gatk_path,
                                   ref = reference_fasta_path,
                                   fai_path = fai_path,
                                   blocksize = blocksize,
                                   cosmic = cosmic_path,
                                   dbsnp = known_snp_vcf_path,
                                   normal_bam_path = normal_bam_path,
                                   tumor_bam_path = tumor_bam_path,
                                   normal_panel = normal_panel,
                                   contEst = contEst,
                                   output_base = os.path.join(step_dir, 'output'))
    )

    outputs = multi_commands(uuid, list(a[0] for a in cmds), thread_count, engine, logger)
    first = True
    with open (out_vcf, "w") as ohandle:
      for cmd, out in cmds:
        with open(out) as handle:
          for line in handle:
            if first or not line.startswith('#'):
              ohandle.write(line)
        first = False
    logger.info('completed running Mutect 2 for: %s and %s' %(normal_bam_path, tumor_bam_path))
