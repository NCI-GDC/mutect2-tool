# mutect2-tool

The GATK3 MuTect2 pipeline employs a "Panel of Normals" to identify additional germline mutations. This panel is generated using TCGA blood normal genomes from thousands of individuals that were curated and confidently assessed to be cancer-free. This method allows for a higher level of confidence to be assigned to somatic variants that were called by the MuTect2 pipeline.

## Installation

```sh
pip install .
```

## Development

* Clone this repository
* Requirements:
  * Python >= 3.9
  * Tox
* `make venv` to create a virtualenv
* `source .venv/bin/activate` to activate new virtualenv
* `make init` to install dependencies and pre-commit hooks
