
contigs="$1"

finishing=AGapEs_analysis/results_AGapEs/finishing/

results=AGapEs_analysis/results_AGapEs/

assembly=AGapEs_analysis/results_AGapEs/assemblies/

map=AGapEs_analysis/results_AGapEs/finishing/ancestral_sequence_map_new


report="AGapEs_analysis/report.txt"

date > $report

echo "----------- Assembly -----------" >> $report
echo "Contig file" >> $report
echo $contigs >> $report
echo "Number of contigs" >> $report
grep ">" $contigs | wc -l  >> $report
echo "total bp in contigs > L" >> $report
grep -v ">" $contigs | tr -d "\n" | wc -c >> $report


echo "" >> $report
echo "--------- Reconstruction --------------" >> $report
echo "" >> $report


echo "number of marker" >> $report
grep ">" ${results}/contigs/families_filtered_uni | wc -l >> $report
echo "number of potential adjacencies" >> $report
wc -l ${results}/scaffold/adjacencies >> $report



echo "Number of CARs" >> $report
grep ">" $map | wc -l >> $report
echo "Length of CARs" >> $report
( grep -A 1 ">" $map | grep -v "CAR" | awk '{print $3}' ; egrep "GAP|MARKER" $map | tail -n 1 | awk '{print $4}' ) | awk 'BEGIN{FS="-"} {print $1}' | grep  "[[:digit:]]" | awk 'p{print $0-p}{p=$0}' >> $report
echo "Total length of CARs" >> $report
( grep -A 1 ">" $map | grep -v "CAR" | awk '{print $3}' ; egrep "GAP|MARKER" $map | tail -n 1 | awk '{print $4}' ) | awk 'BEGIN{FS="-"} {print $1}' | grep  "[[:digit:]]" | awk 'p{print $0-p}{p=$0}' | awk '{sum+=$1} END{print sum}' >> $report
totalCARS=`( grep -A 1 ">" $map | grep -v "CAR" | awk '{print $3}' ; egrep "GAP|MARKER" $map | tail -n 1 | awk '{print $4}' ) | awk 'BEGIN{FS="-"} {print $1}' | grep  "[[:digit:]]" | awk 'p{print $0-p}{p=$0}' | awk '{sum+=$1} END{print sum}'`
echo "Number of gaps per CAR"


echo "Gaps of length 0" >> $report
grep "0-0" ${results}/scaffold/gaps_coordinates_and_length | awk '{print $1}' | cut -c2- > ${results}/gaps_length0
wc -l ${results}/gaps_length0 >> $report



echo "" >> $report
echo "------- Total filled gaps -----------" >> $report
echo "" >> $report





echo "GapFilling Stats consistent" >> $report


echo "Number of consistent filled gaps" >> $report
grep "READS_CONS" ${finishing}/gaps_filled_index_X | wc -l >> $report

echo "" >> $report
echo "Total length of consistent filled gaps" >> $report
grep "READS_CONS" ${finishing}/gaps_filled_index_X | awk '{sum +=$5} END{print sum}' >> $report



echo "" >> $report
echo "GapFilling Stats conflicting" >> $report
echo "Note: contains some of the IS annotated gaps, but we do not count them twice..." >> $report


echo "Number of conflicting filled gaps" >> $report
grep "READS_CONF" ${finishing}/gaps_filled_index_X | wc -l >> $report

echo "Total length of conflicting filled gaps" >> $report
grep "READS_CONF" ${finishing}/gaps_filled_index_X | awk '{sum +=$5} END{print sum}' >> $report



echo "" >> $report
echo "GapFilling Stats IS" >> $report



echo "Number of IS filled gaps " >> $report
wc -l ${results}/assemblies/IS_gaps/gaps_filled >> $report

echo "Total length of IS filled gaps" >> $report
grep "READS_IS" ${finishing}/gaps_filled_index_X | awk '{sum +=$5} END{print sum}' >> $report

echo "" >> $report


echo "Total length of all covered gaps" >> $report
grep "READS" ${finishing}/gaps_filled_index_X | awk '{sum +=$5} END{print sum}' >> $report
cov=`grep "READS" ${finishing}/gaps_filled_index_X | awk '{sum +=$5} END{print sum}'`
echo "Total length of all uncovered gaps" >> $report
grep "FITCH" ${finishing}/gaps_filled_index_X | awk '{sum +=$4} END{print sum}' >> $report



echo "" >> $report
echo "------- Partially filled consistent gaps -----------" >> $report
echo "" >> $report


echo "Uncovered gaps" >> $report
grep "UNCOVERED" ${results}/gapFilling_partial/partGaps_coordinates.index | wc -l >> $report


#dirty hack
cat ${assembly}/consistent/gap_1*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l > outfile
cat ${assembly}/consistent/gap_2*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_3*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_4*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_5*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_6*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_7*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_8*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
cat ${assembly}/consistent/gap_9*.out | grep -v "[[:digit:]]" | egrep -v "A|C|G|T" | wc -l >> outfile
awk '{sum += $1} END {print sum}' outfile >> $report
rm outfile


echo "Partially filled gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial/partGaps_coordinates.index | wc -l >> $report

echo "Total length of partially filled gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial/partGaps_coordinates.index | awk '{sum += $4} END{print sum}' >> $report
echo "Total covered length (coverage > 2) of partially corrected gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial/partGaps_coordinates.index | awk '{sum += $3} END{print sum}' >> $report
part=`grep "PARTIALLY" ${results}/gapFilling_partial/partGaps_coordinates.index | awk '{sum += $3} END{print sum}'`
echo "Total length of uncovered gaps" >> $report
grep "UNCOVERED" ${results}/gapFilling_partial/partGaps_coordinates.index | awk '{sum += $4} END{print sum}' >> $report


echo "" >> $report
echo "------- Partially filled IS gaps -----------" >> $report
echo "" >> $report


echo "Uncovered gaps" >> $report
grep "UNCOVERED" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | wc -l >> $report




echo "Partially filled gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | wc -l >> $report

echo "Total length of partially filled gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | awk '{sum += $4} END{print sum}' >> $report
echo "Total covered length (coverage > 2) of partially corrected gaps" >> $report
grep "PARTIALLY" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | awk '{sum += $3} END{print sum}' >> $report
partIS=`grep "PARTIALLY" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | awk '{sum += $3} END{print sum}'`
echo "Total length of uncovered gaps" >> $report
grep "UNCOVERED" ${results}/gapFilling_partial_IS/partGaps_coordinates.index | awk '{sum += $4} END{print sum}' >> $report




echo "" >> $report

echo "------- Reconstruction -----------" >> $report
echo "" >> $report

echo "Length of reconstructed genome defined by markers" >> $report
grep "MARKER" ${finishing}/ancestral_sequence_map_new | awk '{print $3}' | awk 'BEGIN{FS="-"} {sum += $2-$1+1} END{print sum}' >> $report
mark=`grep "MARKER" ${finishing}/ancestral_sequence_map_new | awk '{print $3}' | awk 'BEGIN{FS="-"} {sum += $2-$1+1} END{print sum}'`
echo "Percentage of reconstructed genome defined by markers" >> $report

echo "scale=8;($mark/$totalCARS)*100" | bc >> $report

echo "Length of reconstructed genome defined by gaps" >> $report
grep "GAP" ${finishing}/ancestral_sequence_map_new | awk '{print $4}' | awk 'BEGIN{FS="-"} {sum += $2-$1+1} END{print sum}' >> $report
gap=`grep "GAP" ${finishing}/ancestral_sequence_map_new | awk '{print $4}' | awk 'BEGIN{FS="-"} {sum += $2-$1+1} END{print sum}'`
echo "scale=8;($gap/$totalCARS)*100" | bc >> $report



echo "Length of reconstructed genome that is supported by reads (= marker + filled gaps + partially filled gaps)" >> $report
echo "$mark+$cov+$part+$partIS" | bc >> $report
echo "Percentage of reconstructed genome that is supported by reads (= marker + filled gaps + partially filled gaps)" >> $report
echo "scale=8;(($mark+$cov+$part+$partIS)/$totalCARS)*100" | bc >> $report











