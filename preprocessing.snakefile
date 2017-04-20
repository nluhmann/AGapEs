configfile: "./config.yaml"

#create directory structure
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/assemblies/consistent')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/finishing/alignments/consistent/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/finishing/alignments/IS/')


rule all:
	input:
		"divide.SUCCESS",
		"IS.SUCCESS"
		




rule fpsac_start:
	input:
		tree=config["tree"],
		contigs=config["contigs"],
		references=config["references"],
	output:			
		"AGapEs_analysis/results_AGapEs/scaffold/adjacencies",
		"AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		"AGapEs_analysis/input/extant_genomes.fasta",
		"AGapEs_analysis/results_AGapEs/contigs/families.fasta",
	log:
		"logs/fpsac/start.log"
	shell:
		"cd FPSAC; ./run_example.sh AGapEs_analysis {input}; cd .."


rule computeGapSequences:
	input:
		ext="AGapEs_analysis/input/extant_genomes.fasta",
		coords="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length"
	output:
		dynamic("AGapEs_analysis/results_AGapEs/finishing/alignments/{gap}.fasta")
	log:
		"logs/computeGapSequences.log"
	shell:
		"cd FPSAC; python2 src/fpsac_extract_extant_gaps_sequences.py ../{input.ext} ../{input.coords} ../AGapEs_analysis/results_AGapEs/finishing/alignments/ ; cd .."


rule prep_IS:
	input:
		gen="AGapEs_analysis/input/extant_genomes.fasta",
		len="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		adjinfo="AGapEs_analysis/adjacencies.info"
	output:
		dynamic("AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_extant")
	log:
		"logs/prepIS.log"
	shell:
		"python2 ./scripts/fpsac_extract_extant_gaps_sequences.py {input.gen} {input.len} AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/ ; python2 ./scripts/extract_IS_extant.py {input.adjinfo} AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/"


rule muscle:
	input: 
		gapf="AGapEs_analysis/results_AGapEs/finishing/alignments/{gap}.fasta",
		tree=config["tree"],
	output:
		"AGapEs_analysis/results_AGapEs/finishing/alignments/{gap}.fasta_ancestral"
	message:
        	"Alignment of {wildcards.gap} to obtain {wildcards.gap}_ancestral"
	log:
		"logs/muscle/{gap}.log"
	shell:
		"cd FPSAC; ./run_muscle.sh {input.gapf} {input.tree}; cd .."




rule muscle2:
	input: 
		gapIS="AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_extant",
		tree=config["tree"],
	output:
		"AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_extant_ancestral"
	message:
        	"Alignment of {wildcards.Igap} to obtain {wildcards.Igap}_ancestral"
	log:
		"logs/muscle2/{Igap}.log"
	shell:
		"cd FPSAC; ./run_muscle.sh {input.gapIS} {input.tree}; cd .."








#todo: generalize from yersinia pestis
rule infos:
	input:
		adj="AGapEs_analysis/results_AGapEs/scaffold/adjacencies",
		iscoords=config["is_coordinates"],
		coords="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
	output:
		"AGapEs_analysis/adjacencies.info",
		"AGapEs_analysis/conflicts.info",
		"AGapEs_analysis/nonConserved.info",
		"AGapEs_analysis/conserved.info",
		"AGapEs_analysis/ISGaps.info"
	log:
		"logs/comp_infos.log"
	shell:
		"python2 scripts/adjacency_infos.py {input.adj} {input.iscoords} {input.coords} {output}"

rule mvIS:
	input:
		"AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_extant_ancestral"
	output:
		"AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_ancestral"
	shell:
		"mv {input} {output}"

rule check_IS:
	input:
		dynamic("AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{Igap}.fasta_ancestral"),
	output:
		"IS.SUCCESS"
	shell:
		"echo '' > IS.SUCCESS"



rule divideGaps_IS:
	input:
		adj="AGapEs_analysis/adjacencies.info",
		conf="AGapEs_analysis/conflicts.info",
		isgaps="AGapEs_analysis/ISGaps.info",
		all=dynamic("AGapEs_analysis/results_AGapEs/finishing/alignments/{gap}.fasta_ancestral")
	output:
		"divide.SUCCESS"
	log:
		"logs/divideGaps.log"
	shell:
		"./scripts/divideGaps.sh {input.adj} {input.conf} {input.isgaps}"



















