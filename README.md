# AGapEs
Gap filling based on template sequence for aDNA data

## Requirements
* python (2.7)

In addition, if available:
* samtools
* bedtools


## Running AGapEs

The following scripts require a template gap sequence and a read mapping in sam format on this template. The assemblies folder should contain an assembly for each gap in fasta format as a concatenation of marker and template gap sequence, and a corresponding read mapping in sam/bam format. The templates folder contains only the templates.

```
run_gapFilling.sh <assemblies folder> <templates folder>
```

Corresponding templates and mappings are indicated by the same gapID. The script provides the required preprocessing of the mappings and calls the central AGapEs implementation for each gap: 

```
python parseAssemblyMappings_Indels.py <sam file> <template> <assembly> <outfile>
```

The repository also contains scripts to fill gaps partially, a bash script will be provided at some point.


