configfile: "./config.yaml"
from os.path import join


IS, = glob_wildcards(join("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/", '{gap,gap_[^/]+}.fasta_ancestral'))



rule all:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partTemplateCorrection.log",
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partTemplateCorrection.log"



rule copyUncovTemplates:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/consistent/copy2.SUCCESS")
	shell:
		"""
		while read p; do cp AGapEs_analysis/results_AGapEs/finishing/alignments/consistent/${{p}}.fasta_ancestral \
		AGapEs_analysis/results_AGapEs/gapFilling_partial/part_templates/${{p}}.fasta_ancestral; done < {input}
		echo "" > AGapEs_analysis/results_AGapEs/assemblies/consistent/copy2.SUCCESS
		"""



rule copyUncovAssemblies:
	input:
		succ="AGapEs_analysis/results_AGapEs/assemblies/consistent/copy2.SUCCESS",
		list="AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/consistent/copyA.SUCCESS")
	shell:
		"""
		while read p; do cp AGapEs_analysis/results_AGapEs/assemblies/consistent/${{p}}.* AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/; done < {input.list}
		echo "" > AGapEs_analysis/results_AGapEs/assemblies/consistent/copyA.SUCCESS
		"""




rule computeStartStopFile:
	input:
		succ="AGapEs_analysis/results_AGapEs/assemblies/consistent/copyA.SUCCESS",
		list="AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled_start_stop"
	shell:
		"while read p; do grep '>' AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/${{p}}.fasta_assembly | awk '{{print $1,$5}}' | sed s/'-'/' '/| tr -d '>' ;  done  < {input.list} > AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled_start_stop"


rule computeUncovCoordinates:
	input:
		list="AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list",
		succ="AGapEs_analysis/results_AGapEs/assemblies/consistent/copyA.SUCCESS"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/merged_uncovered_coordinates_offset"
	shell:
		"""
		rm AGapEs_analysis/results_AGapEs/gapFilling_partial/uncovered_coordinates
		./scripts/getNotCoveredCoordinates.sh AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/ {input.list} AGapEs_analysis/results_AGapEs/gapFilling_partial/
		python2 ./scripts/mergeUncoveredRegions.py AGapEs_analysis/results_AGapEs/gapFilling_partial/uncovered_coordinates > AGapEs_analysis/results_AGapEs/gapFilling_partial/merged_uncovered_coordinates_offset
		"""

rule computePartiallyCovCoords:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled_start_stop",
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/merged_uncovered_coordinates_offset"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates"
	shell:
		"python2 ./scripts/partialGaps_choords.py {input} {output}"


rule runGapFillingPartially:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates"		
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps.log"
	shell:
		"./scripts/run_gapFilling_partially.sh AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/ AGapEs_analysis/results_AGapEs/gapFilling_partial/part_templates/ {input} AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled/ ./scripts/ > AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps.log"
	

rule correctTemplatesPartially:
	input:
		list="AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list",
		log="AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps.log",
		coords="AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates"
											
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial/partTemplateCorrection.log"
	shell:
		"""
		rm AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates.index
		while read GAP; do python2 ./scripts/correctTemplatesPartially.py AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/${{GAP}}.fasta_assembly {input.coords} AGapEs_analysis/results_AGapEs/gapFilling_partial/part_templates/${{GAP}}.fasta_ancestral AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled/${{GAP}}.fasta_ancestral AGapEs_analysis/results_AGapEs/gapFilling_partial/partGaps_coordinates.index AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled/${{GAP}}_*.out ; done < {input.list} > {output}
		"""









rule computeMarkerGapSeq_IS_part:
	input:
		coords="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		align="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/{ISgap}.fasta_ancestral",
		fams="AGapEs_analysis/results_AGapEs/contigs/families.fasta",
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{ISgap}.fasta_assembly"
	log:
		"logs/markerGapSequences_IS.log"
	shell:
		"python2 ./scripts/computeMarkerGapSeq_IS_Fitch.py {wildcards.ISgap} {input.coords} {input.align} {input.fams} AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/"




rule bwa_index:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly"
	output:
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.amb"),
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.ann"),
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.bwt"),
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.pac"),
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.sa")
	log:
		"logs/bwa_index/{gap.log}"
	shell:
		"bwa0_7_9 index {input} > {log}"


rule bwa_map:
	input:
		ref="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly",
		reads=config["reads"],
		index1="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.amb",
		index2="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.ann",
		index3="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.bwt",
		index4="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.pac",
		index5="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.fasta_assembly.sa",
	output:
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam")
	message:
        	"Mapping to {wildcards.gap}"
	log:
		"logs/bwa_map/{gap}.log"
	threads: 1
	shell:
		"(bwa0_7_9 mem -a -t {threads} {input.ref} {input.reads} | "
		"samtools view -Sb -F 4 - > {output}) 2> {log}"


rule samtools_sort:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam"
	output:
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_sorted")
	message:
        	"Sorting {wildcards.gap} bam"
	shell:
		"samtools sort {input} -o {output}"


rule samtools_rmdup:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_sorted"
	output:
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_sorted_dupl")
	message:
        	"Remove duplicates for {wildcards.gap}"
	shell:
		"samtools rmdup -s {input} {output}"

	
rule noClips:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_sorted_dupl"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_noclips"
	message:
        	"Remove clipped mappings for {wildcards.gap}"
	shell:
		"samtools view -h {input} | awk '$6 !~ /H|S/{{print}}' | samtools view -bS - > {output}"

rule computeCov:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_noclips"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bedGraph"
	shell:
		"bedtools genomecov -ibam {input} -bga > {output}"

rule samtools_view:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.bam_noclips"
	output:
		temp("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{gap}.mapped.sam")
	shell:
		"samtools view {input} > {output}"



rule computeStartStopFile_IS:
	input:
		list="AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gaps_notfilled",
		all=expand("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{ISgap}.mapped.sam", ISgap=IS),
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_filled_start_stop"
	shell:
		"while read p; do grep '>' AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/${{p}}.fasta_assembly | awk '{{print $1,$5}}' | sed s/'-'/' '/| tr -d '>' ;  done  < {input.list} > AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_filled_start_stop"


rule computeUncovCoordinates_IS:
	input:
		list="AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gaps_notfilled",
		all=expand("AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/{ISgap}.mapped.sam", ISgap=IS)
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/merged_uncovered_coordinates_offset"
	shell:
		"""
		rm AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/uncovered_coordinates
		./scripts/getNotCoveredCoordinates.sh AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/ {input.list} AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/
		python2 ./scripts/mergeUncoveredRegions.py AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/uncovered_coordinates > AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/merged_uncovered_coordinates_offset
		"""


rule computePartiallyCovCoords_IS:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_filled_start_stop",
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/merged_uncovered_coordinates_offset"
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates"
	shell:
		"python2 ./scripts/partialGaps_choords.py {input} {output}"


rule runGapFillingPartially_IS:
	input:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates"		
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps.log"
	shell:
		"./scripts/run_gapFilling_partially.sh AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/ AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/ {input} AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_filled/ ./scripts/ > AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps.log"


rule correctTemplatesPartially_IS:
	input:
		list="AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gaps_notfilled",
		log="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps.log",
		coords="AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates"
											
	output:
		"AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partTemplateCorrection.log"
	shell:
		"""
		rm AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates.index
		while read GAP; do python2 ./scripts/correctTemplatesPartially.py AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/${{GAP}}.fasta_assembly {input.coords} AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/${{GAP}}.fasta_ancestral AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_filled/${{GAP}}.fasta_ancestral AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/partGaps_coordinates.index AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_filled/${{GAP}}_*.out ; done < {input.list} > {output}
		"""











