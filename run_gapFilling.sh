#assemblies
ASSPATH="$1";
#alignments
ALIPATH="$2"


for ASSEM in ${ASSPATH}*.fasta_assembly
do

	GAP=`basename $ASSEM | cut -f1 -d"." `
	echo $GAP
	
	samtools sort ${ASSPATH}${GAP}.fasta_assembly_mapped.bam ${ASSPATH}${GAP}_s
	samtools rmdup -s ${ASSPATH}${GAP}_s.bam ${ASSPATH}${GAP}_d.bam
	samtools view -h ${ASSPATH}${GAP}_d.bam | awk '$6 !~ /H|S/{print}' | samtools view -bS - > ${ASSPATH}${GAP}_noclips.bam

	bedtools genomecov -ibam ${ASSPATH}${GAP}_noclips.bam -bga > ${ASSPATH}${GAP}.bedGraph

	samtools view ${ASSPATH}${GAP}_noclips.bam > ${ASSPATH}${GAP}_mapped.sam

	python parseAssemblyMappings_Indels.py ${ASSPATH}/${GAP}_mapped.sam $ALIPATH/${GAP}.fasta_ancestral $ASSEM ${ASSPATH}/${GAP}.out

	rm ${ASSPATH}/${GAP}_mapped.sam	
	rm ${ASSPATH}/${GAP}_s.bam
	rm ${ASSPATH}/${GAP}_d.bam
	

done