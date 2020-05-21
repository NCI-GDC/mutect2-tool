"""
Multithreading MuTect2
@author: Shenglai Li
"""

import os
import sys
import time
import glob
import ctypes
import string
import logging
import argparse
import threading
import subprocess
from signal import SIGKILL
from functools import partial
from concurrent.futures import ThreadPoolExecutor


def setup_logger():
    """
    Sets up the logger.
    """
    logger = logging.getLogger("multi_mutect2")
    logger_format = "[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s"
    logger.setLevel(level=logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(logger_format, datefmt="%Y%m%d %H:%M:%S")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def subprocess_commands_pipe(cmd, logger, shell_var=False, lock=threading.Lock()):
    """run pool commands"""
    libc = ctypes.CDLL("libc.so.6")
    pr_set_pdeathsig = ctypes.c_int(1)

    def child_preexec_set_pdeathsig():
        """
        preexec_fn argument for subprocess.Popen,
        it will send a SIGKILL to the child once the parent exits
        """

        def pcallable():
            return libc.prctl(pr_set_pdeathsig, ctypes.c_ulong(SIGKILL))

        return pcallable

    try:
        output = subprocess.Popen(
            shlex.split(cmd),
            shell=shell_var,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            preexec_fn=child_preexec_set_pdeathsig(),
        )
        output.wait()
        with lock:
            logger.info("Running command: %s", cmd)
    except BaseException as e:
        output.kill()
        with lock:
            logger.error("command failed %s", cmd)
            logger.exception(e)
    finally:
        output_stdout, output_stderr = output.communicate()
        with lock:
            logger.error(output_stdout.decode("UTF-8"))
            logger.error(output_stderr.decode("UTF-8"))

def tpe_submit_commands(cmds, thread_count, logger, shell_var=False):
    """run commands on number of threads"""
    with ThreadPoolExecutor(max_workers=thread_count) as e:
        for cmd in cmds:
            e.submit(
                partial(subprocess_commands_pipe, logger=logger, shell_var=shell_var),
                cmd,
            )


def get_region(intervals):
    """get region from intervals"""
    interval_list = []
    with open(intervals, "r") as fh:
        line = fh.readlines()
        for bed in line:
            blocks = bed.rstrip().rsplit("\t")
            intv = "{}:{}-{}".format(blocks[0], int(blocks[1]) + 1, blocks[2])
            interval_list.append(intv)
    return interval_list


def get_file_size(filename):
    """ Gets file size """
    fstats = os.stat(filename)
    return fstats.st_size


def cmd_template(dct):
    """cmd template"""
    lst = [
        "java",
        "-Djava.io.tmpdir=/tmp/job_tmp_${BLOCK_NUM}",
        "-d64",
        "-jar",
        "-Xmx${JAVA_HEAP}",
        "-XX:+UseSerialGC",
        "${GATK_PATH}",
        "-T",
        "MuTect2",
        "-nct",
        "1",
        "-nt",
        "1",
        "-R",
        "${REF}",
        "-L",
        "${REGION}",
        "-I:tumor",
        "${TUMOR_BAM}",
        "-I:normal",
        "${NORMAL_BAM}",
        "--normal_panel",
        "${PON}",
        "--cosmic",
        "${COSMIC}",
        "--dbsnp",
        "${DBSNP}",
        "--contamination_fraction_to_filter",
        "${CONTAMINATION}",
        "-o",
        "${BLOCK_NUM}.mt2.vcf",
        "--output_mode",
        "EMIT_VARIANTS_ONLY",
        "--disable_auto_index_creation_and_locking_when_reading_rods",
    ]
    if dct["dontUseSoftClippedBases"]:
        lst = lst + ["--dontUseSoftClippedBases"]
    template = string.Template(" ".join(lst))
    for i, interval in enumerate(get_region(dct["interval_bed_path"])):
        cmd = template.substitute(
            dict(
                BLOCK_NUM=i,
                JAVA_HEAP=dct["java_heap"],
                GATK_PATH=dct["gatk_path"],
                REF=dct["reference_path"],
                REGION=interval,
                TUMOR_BAM=dct["tumor_bam"],
                NORMAL_BAM=dct["normal_bam"],
                PON=dct["pon"],
                COSMIC=dct["cosmic"],
                DBSNP=dct["dbsnp"],
                CONTAMINATION=dct["contest"],
            )
        )
        yield cmd


def get_args():
    """
    Loads the parser.
    """
    # Main parser
    parser = argparse.ArgumentParser("Internal multithreading MuTect2 calling.")
    # Required flags.
    parser.add_argument(
        "-j",
        "--java_heap",
        required=True,
        help="Java heap memory."
    )
    parser.add_argument(
        "-f",
        "--reference_path",
        required=True,
        help="Reference path."
    )
    parser.add_argument(
        "-r",
        "--interval_bed_path",
        required=True,
        help="Interval bed file."
    )
    parser.add_argument(
        "-t",
        "--tumor_bam",
        required=True,
        help="Tumor bam file."
    )
    parser.add_argument(
        "-n",
        "--normal_bam",
        required=True,
        help="Normal bam file."
    )
    parser.add_argument(
        "-c",
        "--thread_count",
        type=int,
        required=True,
        help="Number of thread."
    )
    parser.add_argument(
        "-p",
        "--pon",
        required=True,
        help="Panel of normals reference path."
    )
    parser.add_argument(
        "-s",
        "--cosmic",
        required=True,
        help="Cosmic reference path."
    )
    parser.add_argument(
        "-d",
        "--dbsnp",
        required=True,
        help="dbSNP reference path."
    )
    parser.add_argument(
        "-e",
        "--contest",
        required=True,
        help="Contamination estimation value from ContEst.",
    )
    parser.add_argument(
        "-m",
        "--dontUseSoftClippedBases",
        action="store_true",
        help="If specified, it will not analyze soft clipped bases in the reads.",
    )
    return parser.parse_args()


def main(args, logger):
    """main"""
    logger.info("Running GATK3.6 MuTect2")
    kwargs = vars(args)
    kwargs["gatk_path"] = "/opt/GenomeAnalysisTK.jar"

    # Start Queue
    tpe_submit_commands(list(cmd_template(kwargs)), kwargs["thread_count"], logger)

    # Check VCFs
    result_vcfs = glob.glob("*.mt2.vcf")
    assert len(result_vcfs) == len(
        get_region(kwargs["interval_bed_path"])
    ), "Missing output!"
    if any(get_file_size(x) == 0 for x in result_vcfs):
        logger.error("Empty VCF detected!")

    # Merge
    first = True
    with open("multi_mutect2_merged.vcf", "w") as oh:
        for out in result_vcfs:
            with open(out) as fh:
                for line in fh:
                    if first or not line.startswith("#"):
                        oh.write(line)
            first = False


if __name__ == "__main__":
    # CLI Entrypoint.
    start = time.time()
    logger_ = setup_logger()
    logger_.info("-" * 80)
    logger_.info("multi_mutect2.py")
    logger_.info("Program Args: %s", " ".join(sys.argv))
    logger_.info("-" * 80)

    args_ = get_args()

    # Process
    logger_.info(
        "Processing tumor and normal bam files %s, %s",
        os.path.basename(args_.tumor_bam),
        os.path.basename(args_.normal_bam),
    )
    main(args_, logger_)

    # Done
    logger_.info("Finished, took %s seconds", round(time.time() - start, 2))
