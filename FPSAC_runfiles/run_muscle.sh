X="$1"
T="$2"

F="../"$X
TREE="../"$T

#echo "---> Aligning gaps sequences with muscle"
python2 src/fpsac_format_fasta.py $F $F"_extant" $F"_tmp"
muscle -in  $F"_extant" -out $F"_muscle_1" -quiet || muscle -in  $F"_extant" -out $F"_muscle_1" -quiet -maxiters 1 || echo "" > $F"_ancestral"


if [ -f $F"_muscle_1" ];
then
python2 src/fpsac_format_fasta.py \
$F"_muscle_1" \
$F"_muscle" \
$F"_tmp"
rm -f $F"_tmp" $F"_muscle_1"

python2 src/fpsac_compute_ancestral_gap_sequence.py \
	$F"_muscle" \
	$TREE \
	$F"_ancestral" || echo "" > $F"_ancestral"
fi


# Output: $FINISHING_ALG/*.fasta_extant - extant gaps FASTA files
# Output: $FINISHING_ALG/*.fasta_muscle - aligned extant gaps FASTA files
# Output: $FINISHING_ALG/*.fasta_ancestral - ancestral gaps FASTA files


