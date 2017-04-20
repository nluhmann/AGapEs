DIR="$1"
TREE="$2"
CONTIGS="$3"
REFS="$4"

./fpsac_without_optimization.sh \
    ../$DIR \
    ../$TREE \
    ../$CONTIGS \
    ../$REFS \
    100 \
    95 \
    1 \
    AGapEs
