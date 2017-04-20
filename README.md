# AGapEs
Gap filling based on template sequence for ancient DNA (aDNA) data [1]

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

### Input
* aDNA reads in fasta/fastq format
* assembled contigs of aDNA reads in fasta format
* assembled reference sequences in fasta format
* phylogenetic tree in newick format, position of the ancient sample marked with @

### Preprocessing

Markers and gap template sequences can be computed using the [FPSAC](https://github.com/cchauve/FPSAC) [2] pipeline. 

```
install_FPSAC.sh
snakemake --snakefile preprocessing.snakefile -j <N>
```

FPSAC computes the template sequences based on a multiple alignment of the respective extant gap sequences. If gaps get too larger because of inadequate marker coverage, this alignment can fail.

### Filling gaps with AGapEs

```
snakemake --snakefile run_gapFilling.snakefile -j <N>
snakemake --snakefile run_partial_gapFilling.snakefile -j <N>
```

### Finishing

```
snakemake --snakefile finishing.snakefile -j <N>
```




# References
[1] Luhmann, Nina, Daniel Doerr, and Cedric Chauve. "Improved assemblies and comparison of two ancient Yersinia pestis genomes." bioRxiv (2017): 073445.

[2] Rajaraman, Ashok, Eric Tannier, and Cedric Chauve. "FPSAC: fast phylogenetic scaffolding of ancient contigs." Bioinformatics 29.23 (2013): 2987-2994.
