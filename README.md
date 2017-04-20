# AGapEs
Gap filling based on template sequence for ancient DNA (aDNA) data

## Requirements

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)
* python
* bwa
* samtools
* bedtools
* FPSAC


# Input and preprocessing

Markers and gap template sequences can be computed using the [FPSAC](https://github.com/cchauve/FPSAC) pipeline. 


# Filling gaps with AGapEs

The following scripts require a template gap sequence and a read mapping in sam format on this template. The assemblies folder should contain an assembly for each gap in fasta format as a concatenation of marker and template gap sequence, and a corresponding read mapping in sam/bam format. The templates folder contains only the templates.

```
run_gapFilling.sh <assemblies folder> <templates folder>
```

Corresponding templates and mappings are indicated by the same gapID. The script provides the required preprocessing of the mappings and calls the central AGapEs implementation for each gap: 

```
python parseAssemblyMappings_Indels.py <sam file> <template> <assembly> <outfile>
```

The repository also contains scripts to fill gaps partially, a bash script will be provided at some point.

# Finishing
