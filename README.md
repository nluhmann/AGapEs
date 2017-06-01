# AGapEs
Gap filling based on template sequence for ancient DNA (aDNA) data [1]

### Requirements

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)
* python2
* bwa
* samtools
* bedtools
* FPSAC (install by running ``` install_FPSAC.sh ```)

### Quickstart ###

```
install_FPSAC.sh
snakemake --snakefile preprocessing.snakefile -j <N>
#Checkpoint: check templates
snakemake --snakefile run_gapFilling.snakefile -j <N>
snakemake --snakefile run_partial_gapFilling.snakefile -j <N>
#Checkpoint: solve conflicts
snakemake --snakefile finishing.snakefile -j <N>
```



## Pipeline
We split the pipeline into 4 snakefiles that should be run successively. If possible, we suggest to allow running the snakefile with multiple cores (parameter -j), since many jobs in the pipeline can be run in parallel.
Each snakefile ends with a "checkpoint" that can be used to e.g. check template sequences or solve conflicts based on gap filling results. The following paragraphs provide detailed information for each pipeline step and how to run the snakefiles.

### Input
* aDNA reads in fasta/fastq format
* assembled contigs of aDNA reads in fasta format (multi-fasta file)
* assembled reference sequences in fasta format (multi-fasta file)
* phylogenetic tree in newick format, placement of the ancient sample marked with @
* IS annotations in extant genomes, format ``` >Genome_ID:start-stop IS_element ```

Before running the pipeline, the location of these input files has to be specified in the file ``` config.yaml```.

### Preprocessing

Markers and gap template sequences can be computed using the [FPSAC](https://github.com/cchauve/FPSAC) [2] pipeline. 

```
install_FPSAC.sh
snakemake --snakefile preprocessing.snakefile -j <N>
```

FPSAC computes the template sequences based on a multiple alignment of the respective extant gap sequences. 

Important files created:
* ``` $DIR/<X>.info ``` - overview files for simple, conflicting and IS gaps. Contains for each potential adjacency: marker extremities, gap ID, in-group and outgrip occurrences and positions, length differences in the extant gaps, IS annotations, conflicting components
* ``` $DIR/results/contigs/families_with_contig_names ``` - markers computed (FPSAC) 
* ``` $DIR/results/contigs/families.fasta ``` - markers sequences (FPSAC)
* ``` $DIR/results/finishing/alignments/ ``` - template sequences for all gaps and underlying extant gap sequences

**CHECKPOINT 1**:
If gaps get too large because of inadequate marker coverage, the multiple alignment of extant gap sequences can fail. In these cases, no template gap sequence will be used in the subsequent gap filling step. If a multiple alignment has been successfully obtained with other software, the file can simply be copied into the respective folder and will be picked up by the pipeline.



### Filling gaps with AGapEs

The two following snakefiles define the full and partial gap filling steps. Since gaps are filled independently, these should be run with multiple threads if possible (-j parameter).

```
snakemake --snakefile run_gapFilling.snakefile -j <N>
snakemake --snakefile run_partial_gapFilling.snakefile -j <N>
```

Important files/directories created:
* ``` $DIR/results/assemblies/<simple|IS|conflicting>/gap_*assembly ``` - template gap sequence together with its flanking marker sequences, for IS-annotated gaps all template alternatives are given. The header of each assembly file is ``` >gapID   leftMarker rightMarker length  gapStart-gapEnd ```
* ``` $DIR/results/assemblies/<simple|IS|conflicting>/gap_*noclips ``` - mapping BAM files
* ``` $DIR/results/assemblies/<simple|IS|conflicting>/gap_*out ``` - gap filling result. The header of each assembly file is ``` >gapID   leftMarker rightMarker length  gapStart-gapEnd dist_to_template``` followed by the reconstructed gap sequence
* ``` $DIR/results/gapFilling_partial/ ``` - assemblies and partial gap filling results for not completely filled simple gaps
* ``` $DIR/results/gapFilling_partial_IS/ ``` - assemblies and partial gap filling results for not completely filled IS gaps

 




**CHECKPOINT 2**: 
Based on the gap filling results, the conflicting components as summarized in $DIR/conflicts.info have to be solved manually. Before running the finishing phase, the two files ``` $DIR/results/scaffold/adjacencies_kept ``` and ``` $DIR/results/scaffold/adjacencies_discarded ``` have to be adjusted accordingly.

### Finishing

The last step will produce a set of scaffolds from kept adjacencies and gap sequences either based on aDNA reads or completed with template sequence. The sequence map also indicates the type of gap and support by aDNA reads.

```
snakemake --snakefile finishing.snakefile -j <N>
```
Important files created:
* ``` $DIR/results/finishing/ancestral_sequence.fasta ``` - resulting ancestral scaffolds
* ``` $DIR/results/finishing/ancestral_sequence_map_augmented ``` - detailed map over the resulting scaffolds, including the status of all gaps included in the assembly. For each gap we have a line:
  	   ``` GAP gapID adjacency position orientation length length_covered_by_reads status OK/DUB ```
We mark all gaps as dubious (DUB) if they are not covered by any reads (COV) or their extant gaps are very unconserved (EXT_LENGTH_DIFF).




## Example

The pipeline has been used to obtain improved assemblies of ancient Yersinia pestis genomes as described in [1]. 
The two config files ```config_London.yaml``` and ```config_Marseille.yaml``` describe all input data for both assemblies. To run the pipeline, download input files from [here](http://paleogenomics.irmacs.sfu.ca/DOWNLOADS/AGAPES_data_results.zip) and the aDNA read files from the respective databases, IDs given in [1].



***


# References
[1] Nina Luhmann, Daniel Doerr, and Cedric Chauve. "Improved assemblies and comparison of two ancient Yersinia pestis genomes." bioRxiv (2017): 073445.

[2] Ashok Rajaraman, Eric Tannier, and Cedric Chauve. "FPSAC: fast phylogenetic scaffolding of ancient contigs." Bioinformatics 29.23 (2013): 2987-2994.
