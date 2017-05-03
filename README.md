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

FPSAC computes the template sequences based on a multiple alignment of the respective extant gap sequences. 

Important files created:
* $DIR/X.info - overview files for simple, conflicting and IS gaps summarizing ingroup and outgroup occurrences, IS annotations and potential conflicting components respectively
* $DIR/results/contigs/families_with_contig_names - markers computed (FPSAC) 
* $DIR/results/contigs/families.fasta - markers fasta sequencs (FPSAC)
* $DIR/results/finishing/alignments/ - template sequences for all gaps and underlying extant gap sequences

CHECKPOINT 1:
If gaps get too large because of inadequate marker coverage, the multiple alignment of extant gap sequences can fail.

### Filling gaps with AGapEs

```
snakemake --snakefile run_gapFilling.snakefile -j <N>
snakemake --snakefile run_partial_gapFilling.snakefile -j <N>
```
Important files/directories created:
* $DIR/results/assemblies/ - assemblies and full gap filling results for all gaps
* $DIR/results/gapFilling_partial/ - assemblies and partial gap filling results for not completely filled simple gaps
* $DIR/results/gapFilling_partial_IS/ - assemblies and partial gap filling results for not completely filled IS gaps

CHECKPOINT 2: 
Based on the gap filling results, the conflicting components have to be solved manually. 

### Finishing

```
snakemake --snakefile finishing.snakefile -j <N>
```
Important files created:
* $DIR/results/finishing/ancestral_sequence_map_augmented
* $DIR/results/finishing/ancestral_sequence.fasta

The last step will produce a set of scaffolds from kept adjacencies and gap sequences either based on aDNA reads or completed with template sequence. We extend the sequence map provided by FPSAC with



# References
[1] Nina Luhmann, Daniel Doerr, and Cedric Chauve. "Improved assemblies and comparison of two ancient Yersinia pestis genomes." bioRxiv (2017): 073445.

[2] Ashok Rajaraman, Eric Tannier, and Cedric Chauve. "FPSAC: fast phylogenetic scaffolding of ancient contigs." Bioinformatics 29.23 (2013): 2987-2994.
