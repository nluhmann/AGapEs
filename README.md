# AGapEs
Gap filling based on template sequence for aDNA data

## Requirements
The following scripts require a template gap sequence and a read mapping in sam format on this template.
In addition, if available:
* python (2.7)
* samtools
* bedtools


## Running AGapEs
```
run_gapFilling.sh <assemblies folder> <templates folder>
```
The assemblies folder should contain for each template the corresponding marker-gap-marker sequence and a mapping of the reads. The templates folder contains only the templates.
Corresponding templates and mappings are indicated by the same gapID. The script provides the required preprocessing of the mappings and calls the central AGapEs implementation for each gap: 

```
python parseAssemblyMappings_Indels.py <sam file> <template> <assembly> <outfile>
```

The repository also contains scripts to fill gaps partially, a bash script will be provided at some point.


