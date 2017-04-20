FINISHING_ALG="$1"
TREE="$2"
SCRIPTS="$3"

echo "---> Aligning gaps sequences with muscle"
GAPS_SEQ=`find ${FINISHING_ALG} -name '*fasta_extant'`
for F in $GAPS_SEQ; do
	if [[ -s $F && ! -s $F"_muscle" ]] ; then
		echo $F
		python $SCRIPTS/fpsac_format_fasta.py \
		$F \
		$F"_extant_form" \
		${FINISHING_ALG}/TMP
		$SCRIPTS/muscle \
		-in  $F"_extant_form" \
		-out $F"_muscle_1" \
		-quiet
		# Below is to handle the well documented muscle bug: http://www.drive5.com/muscle/manual/msa_getletter_bug.html
		if [ ! -f $F"_muscle_1" ];
		then
		echo "-----> Retrying with option maxiters 1: file " $F
		$SCRIPTS/muscle \
		-in  $F"_extant_form" \
		-out $F"_muscle_1" \
		-quiet -maxiters 1
		fi
		python $SCRIPTS/fpsac_format_fasta.py \
		$F"_muscle_1" \
		$F"_muscle" \
		${FINISHING_ALG}/TMP
		rm -f ${FINISHING_ALG}/TMP $F"_muscle_1"
	fi
done
# Output: $FINISHING_ALG/*.fasta_extant - extant gaps FASTA files
# Output: $FINISHING_ALG/*.fasta_muscle - aligned extant gaps FASTA files

# -----------------------------------------------------------
echo
echo "STAGE C.3: COMPUTING ANCESTRAL GAPS SEQUENCES"

echo "---> Computing ancestral gaps sequences"
for F in $GAPS_SEQ; do
	if [[ -s $F ]] ; then
    	python $SCRIPTS/fpsac_compute_ancestral_gap_sequence.py \
		$F"_muscle" \
		$TREE \
		$F"_ancestral"
	fi
done
# Output: $FINISHING_ALG/*.fasta_ancestral - ancestral gaps FASTA files
