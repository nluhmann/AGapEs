# AGapEs
Gap filling based on template sequence for ancient DNA (aDNA) data

### Requirements

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)
* python
* bwa
* samtools
* bedtools
* FPSAC

## Pipeline
We split the pipeline into 4 snakefiles that should be run successively. If possible, we suggest to allow running the snakefile with multiple cores (parameter -j), since many jobs in the pipeline can be run in parallel.
Each snakefile ends with a "checkpoint" that can be used to e.g. check template sequences or solve conflicts based on gap filling results.

## Input
* aDNA reads in fasta/fastq format
* assembled contigs of aDNA reads in fasta format
* assembled reference sequences in fasta format
* phylogenetic tree in newick format, position of the ancient sample marked with @

## Preprocessing

Markers and gap template sequences can be computed using the [FPSAC](https://github.com/cchauve/FPSAC) pipeline. 

```
install_FPSAC.sh
snakemake --snakefile preprocessing.snakefile -j <N>
```

FPSAC computes the template sequences based on a multiple alignment of the respective extant gap sequences. If gaps get too larger because of inadequate marker coverage, this alignment can fail.

## Filling gaps with AGapEs

```
snakemake --snakefile run_gapFilling.snakefile -j <N>
snakemake --snakefile run_partial_gapFilling.snakefile -j <N>
```

## Finishing

```
snakemake --snakefile finishing.snakefile -j <N>
```




# References
