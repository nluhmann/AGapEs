DIR="$1"


echo "Number of IS gaps in directory:"
ls ${DIR}/*assembly | awk 'BEGIN{FS="-"} {print $1'} | sort | uniq | wc -l

echo "Number of gaps in ISInfo:"
grep ">" ../../../ISGaps.info | wc -l

echo "----"

echo "Gaps where all extant gaps are annotated with IS:"
ls ${DIR}/gap_*assembly | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "1 " | awk 'BEGIN{FS="/"} {print $2}'  > gaps_oneVersion
ls ${DIR}/gap_*assembly | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "1 " | awk 'BEGIN{FS="/"} {print $2}' | wc -l

echo "Thereof number of filled gaps:"
while read p; do ls ${DIR}/${p}-*out* | grep "out$" | awk 'BEGIN{FS="-"} {print $1}' | awk 'BEGIN{FS="/"} {print $2}'; done < gaps_oneVersion > gaps_oneVersion_covered
wc -l gaps_oneVersion_covered

echo "CHECK: Thereof number of not filled gaps:"
while read p; do ls ${DIR}/${p}-*out* | grep "uncov$" | awk 'BEGIN{FS="-"} {print $1}' | awk 'BEGIN{FS="/"} {print $2}'; done < gaps_oneVersion > gaps_oneVersion_uncovered
wc -l gaps_oneVersion_uncovered

echo "-----"

echo "Number of gaps with two template versions:"
ls ${DIR}/gap_*assembly | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "2 " | awk 'BEGIN{FS="/"} {print $2}'  > gaps_twoVersion
ls ${DIR}/gap_*assembly | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "2 " | awk 'BEGIN{FS="/"} {print $2}' | wc -l



echo "Gaps where both template versions are covered:"
ls ${DIR}/gap_*out | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "2 " | awk 'BEGIN{FS="-"} {print $1}' | awk 'BEGIN{FS="/"} {print $2}' > gaps_bothVersionCovered
wc -l gaps_bothVersionCovered

echo 'Thereof: IS only annotated in a single extant genome -> take non-IS template for gap'

while read p; do grep ">" ../../finishing/all_extant_IS/${p}-2.fasta_extant | awk '{counter += 1} END{if (counter == 1) print $3}'; done < gaps_bothVersionCovered > gaps_bothVersionCovered_onlyOneExtant
wc -l gaps_bothVersionCovered_onlyOneExtant

echo "CHECK: IS annotated in more than one extant genome"
while read p; do grep ">" ../../finishing/all_extant_IS/${p}-2.fasta_extant | awk '{counter += 1} END{if (counter != 1) print $3}'; done < gaps_bothVersionCovered > gaps_bothVersionCovered_severalExtant
wc -l gaps_bothVersionCovered_severalExtant

echo "Gaps where both template versions are uncovered:"
ls ${DIR}/gap_*uncov | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "2 " | awk 'BEGIN{FS="-"} {print $1}' | awk 'BEGIN{FS="/"} {print $2}' > gaps_bothVersionUncovered
wc -l gaps_bothVersionUncovered



echo "Gaps where only one template version is covered:"
while read p; do ls ${DIR}/${p}-*.out* | grep "out$" | awk 'BEGIN{FS="-"} {print $1}' | sort | uniq -c | grep "1 " | awk 'BEGIN{FS="/"} {print $2}' ; done < gaps_twoVersion > gaps_oneVersionCovered
wc -l gaps_oneVersionCovered


echo "------"

echo "Gaps we consider filled:"
cat gaps_oneVersion_covered gaps_bothVersionCovered_onlyOneExtant gaps_oneVersionCovered > gaps_filled
wc -l gaps_filled


echo 'Gaps not filled, used for partial gapfilling:'
cat gaps_oneVersion_uncovered gaps_bothVersionUncovered gaps_bothVersionCovered_severalExtant > gaps_notfilled
wc -l gaps_notfilled


