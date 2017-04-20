conflicts="$1"
assemblies="$2"

CONFGAPS=`grep "gap_" $conflicts`


for gap in $CONFGAPS
do
	
	if [ -s ${assemblies}/consistent/${gap}.fasta_assembly ];
	then
		mv ${assemblies}/consistent/${gap}.* ${assemblies}/conflicts/
	fi
	if [ -s ${assemblies}/IS_gaps/${gap}-1.fasta_assembly ];
	then
		mv ${assemblies}/IS_gaps/${gap}-* ${assemblies}/conflicts/
	elif [ -s ${assemblies}/IS_gaps/${gap}-2.fasta_assembly ];
	then
		mv ${assemblies}/IS_gaps/${gap}-* ${assemblies}/conflicts/
	fi
done


echo "" > ${assemblies}/conflicts/copy.SUCCESS
