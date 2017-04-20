#assemblies
ASSPATH="$1"
#UncoveredGapsIDs
NOTFILLED="$2"
FOLD="$3"


while read G; do bedtools genomecov -ibam ${ASSPATH}${G}[-.]*bam_noclips -bga > ${ASSPATH}${G}.bedGraph; egrep -w "0$|1$" ${ASSPATH}${G}.bedGraph >> ${FOLD}/uncovered_coordinates; done < $NOTFILLED


