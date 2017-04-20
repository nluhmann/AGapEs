

configfile: "./config.yaml"
from os.path import join

CONS, = glob_wildcards(join("AGapEs_analysis/results_AGapEs/finishing/alignments/consistent/", '{gap,gap_[^/]+}.fasta_ancestral'))
IS, = glob_wildcards(join("AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/", '{gap,gap_[^/]+}.fasta_ancestral'))

snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/assemblies/conflicts/')

snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial/part_templates/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial/part_assemblies/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial/part_filled/')

snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_assemblies/')
snakemake.utils.makedirs('AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_filled/')



rule all:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list",
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/uncov_list",
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/copy.SUCCESS",
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/uncov_list"




rule computeMarkerGapSequences:
	input:
		coords="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		align="AGapEs_analysis/results_AGapEs/finishing/alignments/consistent/{gap}.fasta_ancestral",
		fams="AGapEs_analysis/results_AGapEs/contigs/families.fasta",
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/{gap}.fasta_assembly"
	log:
		"logs/markerGapSequences.log"
	shell:
		"python2 ./scripts/computeMarkerGapSequences.py {wildcards.gap} {input.coords} {input.align} {input.fams} AGapEs_analysis/results_AGapEs/assemblies/consistent/"


rule computeMarkerGapSeq_IS_Fitch:
	input:
		coords="AGapEs_analysis/results_AGapEs/scaffold/gaps_coordinates_and_length",
		align="AGapEs_analysis/results_AGapEs/finishing/all_extant_IS/{ISgap}.fasta_ancestral",
		fams="AGapEs_analysis/results_AGapEs/contigs/families.fasta",
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/{ISgap}.fasta_assembly"
	log:
		"logs/markerGapSequences_IS.log"
	shell:
		"python2 ./scripts/computeMarkerGapSeq_IS_Fitch.py {wildcards.ISgap} {input.coords} {input.align} {input.fams} AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/"



rule bwa_index:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.amb"),
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.ann"),
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.bwt"),
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.pac"),
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.sa")
	log:
		"logs/bwa_index/{gap.log}"
	shell:
		"bwa0_7_9 index {input} > {log}"


rule bwa_map:
	input:
		ref="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly",
		reads=config["reads"],
		index1="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.amb",
		index2="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.ann",
		index3="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.bwt",
		index4="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.pac",
		index5="AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly.sa",
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam")
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
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_sorted")
	message:
        	"Sorting {wildcards.gap} bam"
	shell:
		"samtools sort {input} -o {output}"


rule samtools_rmdup:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_sorted"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_sorted_dupl")
	message:
        	"Remove duplicates for {wildcards.gap}"
	shell:
		"samtools rmdup -s {input} {output}"

	
rule noClips:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_sorted_dupl"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_noclips"
	message:
        	"Remove clipped mappings for {wildcards.gap}"
	shell:
		"samtools view -h {input} | awk '$6 !~ /H|S/{{print}}' | samtools view -bS - > {output}"

rule computeCov:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_noclips"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bedGraph"
	shell:
		"bedtools genomecov -ibam {input} -bga > {output}"

rule samtools_view:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.bam_noclips"
	output:
		temp("AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.mapped.sam")
	shell:
		"samtools view {input} > {output}"


rule gapFilling:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.mapped.sam",
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.fasta_assembly"		
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/{dir}/{gap}.out"
	message:
        	"Run AGapEs for {wildcards.gap}"
	shell:
		"python2 scripts/parseAssemblyMappings_Indels.py {input} {output}"
   

rule moveConflicts:
	input:
		expand("AGapEs_analysis/results_AGapEs/assemblies/consistent/{CONgap}.out", CONgap=CONS),
		expand("AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/{ISgap}.out", ISgap=IS),
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/copy.SUCCESS"
	shell:
		"./scripts/moveConflicts.sh AGapEs_analysis/conflicts.info AGapEs_analysis/results_AGapEs/assemblies/"
		
		


rule moveUncovGaps:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/copy.SUCCESS"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list"
	shell:
		"""
		ls AGapEs_analysis/results_AGapEs/assemblies/consistent/*.out > AGapEs_analysis/results_AGapEs/assemblies/consistent/gap_list
		while read p; do if grep -q "UNCOV" $p; then mv $p $p.uncov; fi ; done < AGapEs_analysis/results_AGapEs/assemblies/consistent/gap_list
		ls AGapEs_analysis/results_AGapEs/assemblies/consistent/*uncov | cut -d "." -f1 | rev | cut -d "/" -f1 | rev > AGapEs_analysis/results_AGapEs/assemblies/consistent/uncov_list
		"""

		
	

rule moveUncovGaps_IS:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/copy.SUCCESS"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/uncov_list"
	shell:
		"""
		ls AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/*.out > AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gap_list
		while read p; do if grep -q "UNCOV" $p; then mv $p $p.uncov; fi ; done < AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gap_list
		ls AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/*uncov | cut -d "." -f1 | rev | cut -d "/" -f1 | rev > AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/uncov_list
		"""

rule moveUncovGaps_conflicts:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/copy.SUCCESS"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/conflicts/uncov_list"
	shell:
		"""
		ls AGapEs_analysis/results_AGapEs/assemblies/conflicts/*.out > AGapEs_analysis/results_AGapEs/assemblies/conflicts/gap_list
		while read p; do if grep -q "UNCOV" $p; then mv $p $p.uncov; fi ; done < AGapEs_analysis/results_AGapEs/assemblies/conflicts/gap_list
		ls AGapEs_analysis/results_AGapEs/assemblies/conflicts/*uncov | cut -d "." -f1 | rev | cut -d "/" -f1 | rev > AGapEs_analysis/results_AGapEs/assemblies/conflicts/uncov_list
		"""


rule ISStats_noclips:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/uncov_list"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gaps_notfilled"
	shell:
		"""
		cd AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/
		../../../../scripts/ISStats_noclips.sh .
		cd ../../../../
		"""


rule copyUncovTemplates_IS:
	input:
		"AGapEs_analysis/results_AGapEs/assemblies/IS_gaps/gaps_notfilled"
	output:
		"AGapEs_analysis/results_AGapEs/assemblies/consistent/copy.SUCCESS"
	shell:
		"""
		while read p; do cp AGapEs_analysis/results_AGapEs/finishing/alignments/IS/${{p}}.fasta_ancestral \
		AGapEs_analysis/results_AGapEs/gapFilling_partial_IS/part_IS_templates/${{p}}.fasta_ancestral; done < {input}
		echo "" > AGapEs_analysis/results_AGapEs/assemblies/consistent/copy.SUCCESS
		"""





