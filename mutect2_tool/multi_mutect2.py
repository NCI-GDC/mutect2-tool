#!/usr/bin/env python3
"""
Multithreading MuTect2
@author: Shenglai Li
"""

import argparse
import concurrent.futures
import logging
import os
import pathlib
import shlex
import subprocess
import sys
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from textwrap import dedent
from types import SimpleNamespace
from typing import IO, Any, Callable, Generator, List, NamedTuple, Optional

logger = logging.getLogger(__name__)

DI = SimpleNamespace(subprocess=subprocess, futures=concurrent.futures,)


class PopenReturn(NamedTuple):
    stderr: str
    stdout: str


CMD_STR = dedent(
    """
    java
    -Djava.io.tmpdir=/tmp/job_tmp_{region_num}
    -d64
    -jar
    -Xmx{java_heap}
    -XX:+UseSerialGC
    {gatk_jar}
    -T MuTect2
    -nct 1
    -nt 1
    -R {reference_path}
    -L {region}
    -I:tumor {tumor_bam}
    -I:normal {normal_bam}
    --normal_panel {pon}
    --cosmic {cosmic}
    --dbsnp {dbsnp}
    --contamination_fraction_to_filter {contamination}
    -o {region_num}.mt2.vcf
    --output_mode
    EMIT_VARIANTS_ONLY
    --disable_auto_index_creation_and_locking_when_reading_rods
    {clipped_bases}
    """
).strip()


def setup_logger():
    """
    Sets up the logger.
    """
    logger_format = "[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s"
    logger.setLevel(level=logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(logger_format, datefmt="%Y%m%d %H:%M:%S")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def subprocess_commands_pipe(cmd: str, timeout: int = 3600, di=DI) -> PopenReturn:
    """run pool commands"""

    output = di.subprocess.Popen(
        shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    try:
        output_stdout, output_stderr = output.communicate(timeout=timeout)
    except Exception:
        output.kill()
        _, output_stderr = output.communicate()
        raise ValueError(output_stderr.decode())

    return PopenReturn(stdout=output_stdout.decode(), stderr=output_stderr.decode(),)


def tpe_submit_commands(
    cmds: List[Any], thread_count: int, fn: Callable = subprocess_commands_pipe, di=DI,
):
    """Run commands on multiple threads.

    Stdout and stderr are logged on function success.
    Exception logged on function failure.
    Accepts:
        cmds (List[str]): List of inputs to pass to each thread.
        thread_count (int): Threads to run
        fn (Callable): Function to run using threads, must accept each element of cmds
    Returns:
        None
    Raises:
        None
    """
    with di.futures.ThreadPoolExecutor(max_workers=thread_count) as executor:
        futures = [executor.submit(fn, cmd) for cmd in cmds]
        for future in di.futures.as_completed(futures):
            try:
                result = future.result()
                logger.info(result.stdout.decode())
                logger.info(result.stderr.decode())
            except Exception as e:
                logger.exception(e)


def yield_bed_regions(intervals_file: str) -> Generator[str, None, None]:
    """Yield region string from BED file."""
    with open(intervals_file, "r") as fh:
        for line in fh:
            chrom, start, end, *_ = line.strip().split()
            interval = "{}:{}-{}".format(chrom, int(start) + 1, end)
            yield interval


def get_file_size(filename: pathlib.Path) -> int:
    """ Gets file size """
    return filename.stat().st_size


def yield_formatted_commands(
    contest: str,
    cosmic: str,
    dbsnp: str,
    gatk_jar: str,
    interval_bed_path: str,
    java_heap: int,
    normal_bam: str,
    not_clipped_bases: bool,
    pon: str,
    reference_path: str,
    tumor_bam: str,
) -> Generator[str, None, None]:
    """Yield commands for each BED interval."""
    clipped_bases = "--dontUseSoftClippedBases" if not_clipped_bases else ""
    for i, region in enumerate(yield_bed_regions(interval_bed_path)):
        cmd = CMD_STR.format(
            region_num=i,
            java_heap=java_heap,
            gatk_jar=gatk_jar,
            reference_path=reference_path,
            region=region,
            tumor_bam=tumor_bam,
            normal_bam=normal_bam,
            pon=pon,
            cosmic=cosmic,
            dbsnp=dbsnp,
            contamination=contest,
            clipped_bases=clipped_bases,
        )
        yield cmd


def setup_parser():
    """
    Loads the parser.
    """
    # Main parser
    parser = argparse.ArgumentParser("Internal multithreading MuTect2 calling.")
    # Required flags.
    parser.add_argument("--java-heap", required=True, help="Java heap memory.")
    parser.add_argument("--reference-path", required=True, help="Reference path.")
    parser.add_argument("--interval-bed-path", required=True, help="Interval bed file.")
    parser.add_argument("--tumor-bam", required=True, help="Tumor bam file.")
    parser.add_argument("--normal-bam", required=True, help="Normal bam file.")
    parser.add_argument(
        "--thread-count", type=int, required=True, help="Number of thread."
    )
    parser.add_argument("--pon", required=True, help="Panel of normals reference path.")
    parser.add_argument("--cosmic", required=True, help="Cosmic reference path.")
    parser.add_argument("--dbsnp", required=True, help="dbSNP reference path.")
    parser.add_argument(
        "--contest", required=True, help="Contamination estimation value from ContEst.",
    )
    parser.add_argument(
        "-m",
        "--dontUseSoftClippedBases",
        action="store_true",
        help="If specified, it will not analyze soft clipped bases in the reads.",
    )
    parser.add_argument(
        "--gatk-jar", default="/usr/local/bin/gatk.jar", required=False,
    )
    return parser


def process_argv(argv: Optional[List] = None) -> namedtuple:
    """Processes argv into namedtuple."""

    parser = setup_parser()

    if argv:
        args, unknown_args = parser.parse_known_args(argv)
    else:
        args, unknown_args = parser.parse_known_args()

    args_dict = vars(args)
    args_dict['extras'] = unknown_args
    run_args = namedtuple('RunArgs', list(args_dict.keys()))
    return run_args(**args_dict)


def merge_files(mutect2_outputs: List[pathlib.Path], out_fh: IO):
    """Write contents of outputs to given file handler."""
    # Merge
    first = True
    for file in mutect2_outputs:
        if get_file_size(file) == 0:
            logger.error("Empty output: %s", file.name)
            continue
        with file.open() as fh:
            for line in fh:
                if first or not line.startswith("#"):
                    out_fh.write(line)
        first = False


def run(run_args):
    """Main script logic.
    Creates muse commands for each BED region and executes in multiple threads.
    """

    run_commands = list(
        yield_formatted_commands(
            run_args.contest,
            run_args.cosmic,
            run_args.dbsnp,
            run_args.gatk_jar,
            run_args.interval_bed_path,
            run_args.java_heap,
            run_args.normal_bam,
            run_args.not_clipped_bases,
            run_args.pon,
            run_args.reference_path,
            run_args.tumor_bam,
        )
    )
    # Start Queue
    tpe_submit_commands(
        run_commands, run_args.thread_count,
    )

    # Check and merge outputs
    p = pathlib.Path('.')
    outputs = list(p.glob("*.mt2.vcf"))

    merged_output_path = "multi_mutect2_merged.vcf"
    with open(merged_output_path, 'w') as fh:
        merge_files(outputs, fh)

    return


def main(argv=None) -> int:
    exit_code = 0

    argv = argv or sys.argv
    args = process_argv(argv)
    setup_logger()

    try:
        run(args)
    except Exception as e:
        logger.exception(e)
        exit_code = 1
    return exit_code


if __name__ == "__main__":
    # CLI Entrypoint.
    retcode = 0

    try:
        retcode = main()
    except Exception as e:
        retcode = 1
        logger.exception(e)

    sys.exit(retcode)


# __END__
