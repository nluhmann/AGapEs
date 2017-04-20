configfile: "./config.yaml"
from os.path import join


snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/finishing/alignments_final/')


#Todo 1: solve conflicts




rule all:
	input:
		"AGapEs_analysis/report.txt",
		"AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new_augmented"



rule sortAndIndex:
	input:
		"AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		"AGapEs_analysis/results_AGapEs/finishing/",
		"AGapEs_analysis/results_AGapEs/assemblies/"
	output:
		"AGapEs_analysis/results_AGapEs/finishing/gaps_filled_index_X"
	shell:
		"""
		python2 ./scripts/index_ancestral_gap_sequences.py {input} X 
		python2 ./scripts/sort_ancestral_gap_sequences.py {input}
		cp AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled/*ancestral AGapEs_analysis/results_AGapEs/finishing/alignments_final/
		cp AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_filled/*ancestral AGapEs_analysis/results_AGapEs/finishing/alignments_final/
		"""
		




rule run_FPSAC_finishing:
	input:
		tree=config["tree"],
		contigs=config["contigs"],
		references=config["references"],
		index="AGapEs_analysis/results_AGapEs/finishing/gaps_filled_index_X",
	output:
		"AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new"
	log:
		"logs/fpsac/finishing.log"
	shell:
		"cd FPSAC; ./run_finishing.sh AGapEs_analysis {input.tree} {input.contigs} {input.references}; cd .."



rule extendAncestralMap:
	input:
		"AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new",
		"AGapEs_analysis/results_AGapEs/finishing/gaps_filled_index_X",
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates.index",
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates.index",
		"AGapEs_analysis/adjacencies.info"
	output:
		"AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new_augmented"
	shell:
		"python2 ./scripts/extending_ancestral_map.py {input}"


rule createReport:
	input:
		contigs=config["contigs"],
		map="AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new"
	output:
		"AGapEs_analysis/report.txt"
	shell:
		"./scripts/createReport.sh {input.contigs}"






