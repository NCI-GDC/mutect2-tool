# GDC GATK3 MuTect2
![Version badge](https://img.shields.io/badge/GATK3.6-nightly--2016--02--25--gf39d340-<COLOR>.svg)

The GATK3 MuTect2 pipeline employs a "Panel of Normals" to identify additional germline mutations. This panel is generated using TCGA blood normal genomes from thousands of individuals that were curated and confidently assessed to be cancer-free. This method allows for a higher level of confidence to be assigned to somatic variants that were called by the MuTect2 pipeline.

Original GATK3 MuTect2: https://gatkforums.broadinstitute.org/gatk/discussion/9183/how-to-call-somatic-snvs-and-indels-using-mutect2

## GATK3

Important note:

* The GDC GATK MuTect2 version was frozen to the version when we delivered our first data release. GATK team normally do not keep nightly version beyond 30 days, so that it makes really difficult to re-build the identical docker image.<br>
However, according to GATK team, it seems reasonable to use GATK3.7 as a replacement.<br>
https://gatkforums.broadinstitute.org/gatk/discussion/9406/where-can-i-find-the-gdc-mutect2-version
* Please contact GATK team for the GATK3.7 `GenomeAnalysisTK.jar`.

## Docker

There are two `Dockerfile`s for different purposes:

* Vanilla MuTect2
  * `/docker/Dockerfile` : `GATK3.6, nightly-2016-02-25-gf39d340` docker without additional features.
* Multi-threading MuTect2
  * `/docker/multi_mutect2/Dockerfile` : A python multi-threading implementation on MuTect2 function. Achieve `scatter/gather` method on Docker level.

## How to build

https://docs.docker.com/engine/reference/builder/

The docker images are tested under multiple environments. The most tested ones are:
* Docker version 19.03.2, build 6a30dfc
* Docker version 18.09.1, build 4c52b90
* Docker version 18.03.0-ce, build 0520e24
* Docker version 17.12.1-ce, build 7390fc6

## For external users

There is a production-ready CWL example at https://github.com/NCI-GDC/mutect2-cwl which uses the docker images that are built from the `Dockerfile`s in this repo.

To run multi-threading GATK3 MuTect2:

```
[INFO] [20200109 04:10:13] [multi_mutect2] - --------------------------------------------------------------------------------
[INFO] [20200109 04:10:13] [multi_mutect2] - multi_mutect2.py
[INFO] [20200109 04:10:13] [multi_mutect2] - Program Args: docker/multi_mutect2/multi_mutect2.py -h
[INFO] [20200109 04:10:13] [multi_mutect2] - --------------------------------------------------------------------------------
usage: Internal multithreading MuTect2 calling. [-h] -j JAVA_HEAP -f
                                                REFERENCE_PATH -r
                                                INTERVAL_BED_PATH -t TUMOR_BAM
                                                -n NORMAL_BAM -c THREAD_COUNT
                                                -p PON -s COSMIC -d DBSNP -e
                                                CONTEST [-m]

optional arguments:
  -h, --help            show this help message and exit
  -j JAVA_HEAP, --java_heap JAVA_HEAP
                        Java heap memory.
  -f REFERENCE_PATH, --reference_path REFERENCE_PATH
                        Reference path.
  -r INTERVAL_BED_PATH, --interval_bed_path INTERVAL_BED_PATH
                        Interval bed file.
  -t TUMOR_BAM, --tumor_bam TUMOR_BAM
                        Tumor bam file.
  -n NORMAL_BAM, --normal_bam NORMAL_BAM
                        Normal bam file.
  -c THREAD_COUNT, --thread_count THREAD_COUNT
                        Number of thread.
  -p PON, --pon PON     Panel of normals reference path.
  -s COSMIC, --cosmic COSMIC
                        Cosmic reference path.
  -d DBSNP, --dbsnp DBSNP
                        dbSNP reference path.
  -e CONTEST, --contest CONTEST
                        Contamination estimation value from ContEst.
  -m, --dontUseSoftClippedBases
                        If specified, it will not analyze soft clipped bases
                        in the reads.
```

## For GDC users

See https://github.com/NCI-GDC/gdc-somatic-variant-calling-workflow.
