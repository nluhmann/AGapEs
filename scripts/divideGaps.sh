adjacency="$1"
conflicts="$2"
IS="$3"

#mkdir ${alignments}/consistent
#mkdir ${alignments}/conflicts
#mkdir ${alignments}/IS

ISGAPS=`grep "gap_" $IS`
CONFGAPS=`grep "gap_" $conflicts`
ADJ=`grep "gap_" $adjacency`

alignments="AGapEs_analysis/results_AGapEs/finishing/alignments"

#for gap in $CONFGAPS
#do
	
#	if [ -s ${alignments}/${gap}.fasta ];
#	then
#		cp ${alignments}/${gap}.*ancestral ${alignments}/conflicts
#	else
#		echo "no file for ${gap} to move to conflicts"
#	fi
#done



for gap in $ISGAPS
do		
	if [ -s ${alignments}/${gap}.fasta ];
	then
		cp $alignments/${gap}.*ancestral ${alignments}/IS
	else
		echo "no file for ${gap} to move to IS"
	fi
done


for gap in $ADJ
do		
	if [ -s ${alignments}/${gap}.fasta ];
	then
		if [ ! -s ${alignments}/conflicts/${gap}.*ancestral ];
		then
			if [ ! -s ${alignments}/IS/${gap}.*ancestral ];
			then
				cp $alignments/${gap}.*ancestral ${alignments}/consistent
			else
				echo "no file for ${gap} to move to consistent"
			fi
		fi
	fi
done

echo "" > divide.SUCCESS
