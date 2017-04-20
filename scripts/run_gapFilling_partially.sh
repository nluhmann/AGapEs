#assemblies
ASSPATH="$1";
#alignments
ALIPATH="$2";
#partGaps_coordinates
GAPCOORDS="$3";
#outdir
OUT="$4"
#scripts
SCR="$5"

while read LINE
do
	ary=($LINE)
	GAP=${ary[0]}
	echo $GAP
	COOR=${ary[1]}
	echo $COOR
	MODE=${ary[2]}
	samtools index ${ASSPATH}/${GAP}[-.]*bam_noclips
	samtools view ${ASSPATH}/${GAP}[-.]*bam_noclips ${GAP}:${COOR} > ${ASSPATH}/${GAP}_${COOR}
	
	python2 ${SCR}/parseAssemblyMappings_partially.py ${ASSPATH}/${GAP}_${COOR} $ALIPATH/${GAP}.fasta_ancestral ${ASSPATH}/${GAP}[-.]*fasta_assembly ${MODE} ${OUT}/${GAP}_${COOR}.out


done < $GAPCOORDS








# while read GAP
# do
# 	grep ${GAP} ${MERGE} > tmp_offset
# 	MAP=`ls ${ASSPATH}/${GAP}_[0-9]*| sort -n -t "-" -k2`	
# 	ary=($MAP)
# 	LEN=${#ary[@]}
# 	var=$((${LEN}-1))
# 	for key in "${!ary[@]}"; do
# 		echo "$key ${ary[$key]}"
# 		if [ $key -eq 0 ]; then
# 			echo "LEFT"
# 			python /Users/Nina/Desktop/Pestis\ Analysis\ Scripts/parseAssemblyMappings_partially.py ${ary[$key]} $ALIPATH/${GAP}.fasta_ancestral ${ASSPATH}/${GAP}.fasta_assembly "left" ${OUT}/${GAP}_${key}.out
# 		elif [ $key -eq $var ]; then
# 			echo "RIGHT"
# 			python /Users/Nina/Desktop/Pestis\ Analysis\ Scripts/parseAssemblyMappings_partially.py ${ary[$key]} $ALIPATH/${GAP}.fasta_ancestral ${ASSPATH}/${GAP}.fasta_assembly "right" ${OUT}/${GAP}_${key}.out
# 		else
# 			echo "INTERNAL"
# 			python /Users/Nina/Desktop/Pestis\ Analysis\ Scripts/parseAssemblyMappings_partially.py ${ary[$key]} $ALIPATH/${GAP}.fasta_ancestral ${ASSPATH}/${GAP}.fasta_assembly "internal" ${OUT}/${GAP}_${key}.out
# 		fi
# 	done
# 	
# 	rm tmp_offset
# 
# done < $UNCOV
