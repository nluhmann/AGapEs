# Modification 2015: When computing ancestral gap sequence, use all adjacencies instead of kept ones




# FPSAC, Fast Phylogenetic Scaffolding of Ancient Contigs.
# June 2013. 
# cedric.chauve@sfu.ca

# PARAMETERS
#
# $1 = directory where data and results are stored, directory name (created if required)
# $2 = species tree with marked ancestor, file name
#      file format: Nexus file, with no annotation, with branch lengths and ancestor node marked by @
# $3 = ancient contigs sequence file
#      file format: FASTA
# $4 = extant genomes sequence file
#      file format: FASTA
# $5 = length threshold used to define homologous families, integer
# $6 = similarity threshold used to define homologous families, integer
# $7 = upper bound on ancestral multiplicity used in the scaffolding phase, integer
# $8 = ancestor name, single word text

date
uname -a
echo

# -----------------------------------------------------------
# DIRECTORIES AND FILE NAMES
# -----------------------------------------------------------
echo "CREATING DIRECTORY STRUCTURE AND FORMATTING INPUT FILES"

DIR=$1
DATA=$DIR/input
RESULTS=$DIR/results_$8
CONTIGS=$RESULTS/contigs
SCAFFOLD=$RESULTS/scaffold
FINISHING=$RESULTS/finishing
FINISHING_ALG=$FINISHING/alignments_final

LG_THRESHOLD=$5
SIM_THRESHOLD=$6
MULT_THRESHOLD=$7
ANCESTOR_NAME=$8


echo
echo "DIRECTORY: $DIR"
echo "ANCESTOR NAME: $ANCESTOR_NAME"
echo

# echo "---> Creating/cleaning directories"
# mkdir -p $DATA
# mkdir -p $CONTIGS
# mkdir -p $SCAFFOLD
# mkdir -p $FINISHING_ALG
# rm    -f $FINISHING_ALG/*fasta_muscle

TREE=$DATA/tree
SPECIES=$DATA/species
SPECIES_PAIRS=$DATA/species_pairs

EXTANT_GENOMES_FA=$DATA/extant_genomes.fasta 
EXTANT_GENOMES_LG=$DATA/extant_genomes_length
EXTANT_GENOMES_DB=$DATA/extant_genomes_db

ANCESTRAL_CONTIGS_FA=$DATA/ancestral_contigs.fasta 
ANCESTRAL_CONTIGS_LG=$DATA/ancestral_contigs_length

MEGABLAST_HITS=$CONTIGS/megablast_hits
ANCIENT_EXTANT_HITS=$CONTIGS/ancient_extant_hits

FAMILIES_COORDINATES_ALL=$CONTIGS/families_with_contig_names
FAMILIES_PROFILES_ALL=$CONTIGS/families_profiles
FAMILIES_ERRORS=$CONTIGS/families_errors
FAMILIES_ANCESTRAL_CONTENT=$CONTIGS/families_ancestral_content
FAMILIES_COORDINATES_FILTERED=$CONTIGS/families_filtered
FAMILIES_COORDINATES_FILTERED_UNI=$CONTIGS/families_filtered_uni
FAMILIES_COORDINATES_FILTERED_DOUBLED=$CONTIGS/families_filtered_doubled
FAMILIES_FA=$CONTIGS/families.fasta
FAMILIES_LG=$CONTIGS/families_length

ADJACENCIES_ALL=$SCAFFOLD/adjacencies
ADJACENCIES_WGTD=$SCAFFOLD/adjacencies_weighted
ADJACENCIES_KEPT=$SCAFFOLD/adjacencies_kept
ADJACENCIES_DISC=$SCAFFOLD/adjacencies_discarded
ADJACENCIES_CARS_ALL=$SCAFFOLD/adjacencies_CARS
RSI_ALL=$SCAFFOLD/repeat_spanning_intervals
RSI_WGTD=$SCAFFOLD/repeat_spanning_intervals_weighted
RSI_KEPT=$SCAFFOLD/repeat_spanning_intervals_kept
RSI_DISC=$SCAFFOLD/repeat_spanning_intervals_disc

REPEAT_CLUSTERS=$SCAFFOLD/repeat_clusters

SCAFFOLD_ORDER_DOUBLED=$SCAFFOLD/scaffold_order_doubled_new
SCAFFOLD_ORDER=$SCAFFOLD/scaffold_order_new
ADJACENCIES_UNASSIGNED=$SCAFFOLD/adjacencies_unassigned_new
SCAFFOLD_FA=$SCAFFOLD/scaffold_new.fasta
SCAFFOLD_FINISHED_FA=$FINISHING/ancestral_sequence_new.fasta
SCAFFOLD_MAP=$FINISHING/ancestral_sequence_map_new

GAPS_COORDINATES=$SCAFFOLD/gaps_coordinates
GAPS_COORDINATES_LG=$SCAFFOLD/gaps_coordinates_and_length

REPORT=$DIR/report_$8

# echo "---> Formatting and copying input files"
# python src/fpsac_format_fasta.py $4 $EXTANT_GENOMES_FA    $EXTANT_GENOMES_LG
# python src/fpsac_format_fasta.py $3 $ANCESTRAL_CONTIGS_FA $ANCESTRAL_CONTIGS_LG
# cp $2 $TREE
# 
# # -----------------------------------------------------------
# # STAGE A. COMPUTING HOMOLOGOUS FAMILIES
# # -----------------------------------------------------------
# echo
# #-----------------------------------------------------------
# echo "STAGE A.1: COMPUTING CONTIGS/EXTANT GENOMES HITS"
# 
# echo "---> Computing hits: BLAST index for extant genomes (makembindex)"
# makembindex \
#     -input  ${EXTANT_GENOMES_FA} \
#     -output ${EXTANT_GENOMES_DB} \
#     -verbosity quiet
# echo "---> Computing hits: BLAST database for extant genomes (makeblastdb)"
# makeblastdb \
#     -in     ${EXTANT_GENOMES_FA} \
#     -dbtype nucl 
# echo "---> Computing hits (blastn -task megablast)"
# blastn -task megablast \
#     -db         ${EXTANT_GENOMES_FA} \
#     -query      ${ANCESTRAL_CONTIGS_FA} \
#     -use_index  True \
#     -index_name ${EXTANT_GENOMES_DB} \
#     -out        ${MEGABLAST_HITS} \
#     -outfmt 6 
# echo "---> Formatting hits"
# python src/fpsac_format_blast_hits.py \
#     ${MEGABLAST_HITS} \
#     ${ANCESTRAL_CONTIGS_LG} \
#     ${ANCIENT_EXTANT_HITS}
# # Output: $ANCIENT_EXTANT_HITS - contigs/genomes hits file
# # 
# # -----------------------------------------------------------
# # echo
# # echo "STAGE A.2: COMPUTING HOMOLOGOUS FAMILIES AND PROFILES"
# # 
# echo "---> Computing families"
# python src/fpsac_compute_families_coordinates_and_profiles.py \
#     ${ANCIENT_EXTANT_HITS} \
#     ${FAMILIES_COORDINATES_ALL}_1 \
#     ${FAMILIES_PROFILES_ALL}_1 \
#     ${LG_THRESHOLD} \
#     ${SIM_THRESHOLD}
# python src/fpsac_correct_families.py \
#      ${FAMILIES_COORDINATES_ALL}_1 \
#      ${FAMILIES_PROFILES_ALL}_1 \
#      ${FAMILIES_COORDINATES_ALL} \
#      ${FAMILIES_PROFILES_ALL} \
#      ${FAMILIES_ERRORS}
# rm -f ${FAMILIES_COORDINATES_ALL}_1 ${FAMILIES_PROFILES_ALL}_1
# # Output: ${FAMILIES_COORDINATES_ALL} - families and markers coordinates
# # Output: ${FAMILIES_PROFILES_ALL}    - families extant profiles
# # Output: ${FAMILIES_ERRORS}          - families where markers could not be oriented (discarded from further analyses)
# # 
# # -----------------------------------------------------------
# # echo
# # echo "STAGE A.3: COMPUTING ANCESTRAL CONTENT/MULTIPLICITIES"
# # 
# echo "---> Computing ancestral content/multiplicities"
# python src/fpsac_compute_families_ancestral_content.py \
#     $TREE \
#     ${FAMILIES_PROFILES_ALL} \
#     ${FAMILIES_ANCESTRAL_CONTENT}
# # Output: ${FAMILIES_ANCESTRAL_CONTENT} - families ancestral content
# 
# # -----------------------------------------------------------
# # STAGE B. SCAFFOLDING
# # -----------------------------------------------------------
# echo
# echo "STAGE B.1: SPECIES, SPECIES PAIRS, DOUBLING MARKERS, FILTERING HIGH MULTIPLICITIES FAMILIES"
# 
# echo "---> Computing species pairs list"
# python src/fpsac_compute_species_pairs_list.py \
#    $TREE \
#    ${SPECIES_PAIRS}
# # Output: ${SPECIES_PAIRS} - pairs of species whose genomes are compared
# 
# echo "---> Computing species list"
# python src/fpsac_compute_species_list.py \
#     $TREE \
#     ALL \
#     $SPECIES
# # # Output: $SPECIES - list of species in the tree
# # 
#  echo "---> Filtering high multiplicity homologous families"
# python src/fpsac_filter_families.py \
#     ${FAMILIES_COORDINATES_ALL} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     ${MULT_THRESHOLD} \
#     ${FAMILIES_COORDINATES_FILTERED}
# # # Output: ${FAMILIES_COORDINATES_FILTERED} - families with multiplicity below chosen threshold
# 
# 
# #### changed by Nina
#  echo "---> Filtering high multiplicity homologous families"
# python filterUniversalFamilies.py \
# 	${FAMILIES_COORDINATES_FILTERED} \
# 	${FAMILIES_COORDINATES_FILTERED_UNI}
# # # Output: ${FAMILIES_COORDINATES_FILTERED_UNI} - universal families
# 
# 
# 
# echo "---> Doubling markers"
# python src/fpsac_double_markers.py \
#     ${FAMILIES_COORDINATES_FILTERED_UNI} \
#     ${FAMILIES_COORDINATES_FILTERED_DOUBLED}
# # Output: ${FAMILIES_COORDINATES_FILTERED_DOUBLED} - doubled markers
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE B.2.a: COMPUTING AND WEIGHTING ADJACENCIES"
# 
# echo "---> Computing adjacencies"
# python src/fpsac_compute_conserved_syntenies.py \
#     ${FAMILIES_COORDINATES_FILTERED_DOUBLED} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     $SPECIES \
#     ${SPECIES_PAIRS} \
#     1 \
#     ADJ \
#     ${ADJACENCIES_ALL}
# # Note: parameter 1 indicates that the considered genomes are circular
# # Output: $ADJACENCIES_ALL - Dollo conserved adjacencies between markers
# 
# echo "---> Weighting adjacencies"
# python src/fpsac_weight.py \
#     ${ADJACENCIES_ALL} \
#     $TREE \
#     ${ADJACENCIES_WGTD} \
#     d
# # Output: $ADJACENCIES_WGTD - weighted adjacencies
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE B.2.b: COMPUTING AND WEIGHTING REPEAT SPANNING INTERVALS"
# 
# echo "---> Computing repeat spanning intervals"
# python src/fpsac_compute_conserved_syntenies.py \
#     ${FAMILIES_COORDINATES_FILTERED_DOUBLED} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     $SPECIES \
#     ${SPECIES_PAIRS} \
#     1 \
#     RSI \
#     ${RSI_ALL}
# # Output: ${RSI_ALL} - repeat spanning intervals
# 
# echo "---> Weighting repeat spanning intervals"
# python src/fpsac_weight.py \
#     ${RSI_ALL} \
#     $TREE \
#     ${RSI_WGTD}_tmp \
#     d
# python src/fpsac_correct_weighted_repeat_spanning_intervals.py \
#     ${RSI_WGTD}_tmp \
#     ${RSI_ALL} \
#     ${RSI_WGTD}
# rm -f  ${RSI_WGTD}_tmp
# # Output: $RSI_WGTD - weighted repeat spanning intervals
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE B.3.a: LINEARIZATION/CIRCULARIZATION"
# 
# echo "---> Filtering adjacencies"
# python src/fpsac_filter_adjacencies.py \
#     ${ADJACENCIES_WGTD} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     1 \
#     ${ADJACENCIES_KEPT} \
#     ${ADJACENCIES_DISC} 
# # Note:  parameter 1 indicates that the considered genomes are circular
# # Output: $ADJACENCIES_KEPT - kept adjacencies
# # Output: $ADJACENCIES_DISC - discarded adjacencies
# 
# echo "---> Filtering repeat spanning intervals"
# python src/fpsac_filter_repeat_spanning_intervals.py \
#     ${RSI_WGTD} \
#     ${ADJACENCIES_KEPT} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     ${RSI_KEPT} \
#     ${RSI_DISC}
# # Output: $RSI_KEPT - kept repeat spanning intervals
# # Output: $RSI_DISC - discarded repeat spanning intervals
# 
# # -----------------------------------------------------------
# echo
echo "STAGE B.3.b: COMPUTING SCAFFOLDS FROM FILTERED ADJACENCIES AND REPEAT SPANNING INTERVALS"

echo "---> Computing doubled scaffolds"
python2 src/fpsac_compute_scaffold_order.py \
    ${RSI_KEPT} \
    ${ADJACENCIES_KEPT} \
    ${FAMILIES_ANCESTRAL_CONTENT} \
    ${SCAFFOLD_ORDER_DOUBLED} \
    ${ADJACENCIES_UNASSIGNED} \
    ${ANCESTOR_NAME}
# Output: $SCAFFOLD_ORDER_DOUBLED  - scaffold order with doubled markers
# Output: $ADJACENCIES_UNASSIGNED - adjacencies unassigned to any scaffold

echo "---> Undoubling markers"
python2 src/fpsac_halve_scaffold.py \
    ${SCAFFOLD_ORDER_DOUBLED} \
    ${SCAFFOLD_ORDER}
# Output: $SCAFFOLD_ORDER - scaffold order

# -----------------------------------------------------------
# echo
# echo "STAGE B.3.c: COMPUTING ADJACENCIES BETWEEN CARS AND REPEAT CLUSTERS"
# 
# echo "---> Misc: Computing adjacencies between CARs extremities"
# python src/fpsac_compute_cars_adjacencies.py \
#     ${FAMILIES_COORDINATES_FILTERED_DOUBLED} \
#     ${SCAFFOLD_ORDER_DOUBLED} \
#     $SPECIES \
#     1 \
#     ${ADJACENCIES_CARS_ALL}
# # Note: parameter 1 indicates that the considered genomes are circular
# # Output: $ADJACENCIES_CARS_ALL - all adjacencies between CARS extremities
# 
# echo "---> Misc: Computing maximal repeat clusters"
# python src/fpsac_compute_repeat_clusters.py \
#     ${ADJACENCIES_ALL} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     ${REPEAT_CLUSTERS}
# # Output: $REPEAT_CLUSTERS - list of all maximal repeat clusters
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE C.1: COMPUTING ANCESTRAL GAPS LENGTHS"
# 
# echo "---> Extracting conserved extant gaps"
# python src/fpsac_extract_gaps.py \
#     ${ADJACENCIES_ALL} \
#     ${RSI_KEPT} \
#     ${FAMILIES_ANCESTRAL_CONTENT} \
#     ${GAPS_COORDINATES}
# # Output: $GAPS_COORDINATES - coordinates in extant genomes of all gaps present in the scaffolds
# 
# echo "---> Computing gaps length"
# python src/fpsac_compute_supported_gaps.py \
#     ${GAPS_COORDINATES} \
#     ${SPECIES_PAIRS} \
#     1 0 \
#     ${GAPS_COORDINATES_LG}
# # Output: $GAPS_COORDINATES_LG - coordinates in extant genomes and length in the ancestral genome of Dollo conserved gaps
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE C.2: COMPUTING AND ALIGNING EXTANT GAPS SEQUENCES"
# 
# echo "---> Computing gaps sequences"
# python src/fpsac_extract_extant_gaps_sequences.py \
#     ${EXTANT_GENOMES_FA} \
#     ${GAPS_COORDINATES_LG} \
#     ${FINISHING_ALG}/
# 
# echo "---> Aligning gaps sequences with muscle"
# GAPS_SEQ=`find ${FINISHING_ALG} -name '*fasta'`
# for F in $GAPS_SEQ; do
#     python src/fpsac_format_fasta.py \
# 	$F \
# 	$F"_extant" \
# 	$FINISHING/TMP
#     muscle \
# 	-in  $F"_extant" \
# 	-out $F"_muscle_1" \
# 	-quiet
#     # Below is to handle the well documented muscle bug: http://www.drive5.com/muscle/manual/msa_getletter_bug.html
#     if [ ! -f $F"_muscle_1" ];
#     then
# 	echo "-----> Retrying with option maxiters 1: file " $F
# 	muscle \
# 	-in  $F"_extant" \
# 	-out $F"_muscle_1" \
# 	-quiet -maxiters 1
#     fi
#     python src/fpsac_format_fasta.py \
# 	$F"_muscle_1" \
# 	$F"_muscle" \
# 	$FINISHING/TMP
#     rm -f $FINISHING/TMP $F"_muscle_1"
# done
# # Output: $FINISHING_ALG/*.fasta_extant - extant gaps FASTA files
# # Output: $FINISHING_ALG/*.fasta_muscle - aligned extant gaps FASTA files
# 
# # -----------------------------------------------------------
# echo
# echo "STAGE C.3: COMPUTING ANCESTRAL GAPS SEQUENCES"
# 
# echo "---> Computing ancestral gaps sequences"
# for F in $GAPS_SEQ; do
#     python src/fpsac_compute_ancestral_gap_sequence.py \
# 	$F"_muscle" \
# 	$TREE \
# 	$F"_ancestral"
# done
# Output: $FINISHING_ALG/*.fasta_ancestral - ancestral gaps FASTA files

# -----------------------------------------------------------
echo
echo "STAGE C.4: COMPUTING ANCESTRAL SCAFFOLD SEQUENCE AND MAP"

# echo "---> Computing sequences of ancestral markers"
# python src/fpsac_extract_families_sequences.py \
#     ${FAMILIES_COORDINATES_ALL} \
#     ${ANCESTRAL_CONTIGS_FA} \
#     ${FAMILIES_FA} \
#     ${FAMILIES_LG}
# Output: $FAMILIES_FA - FASTA sequences of ancestral markers
# Output: $FAMILIES_LG - length of sequences of ancestral markers

echo "---> Computing scaffold sequence with Ns in gaps"
python2 src/fpsac_compute_scaffold_sequence_N.py \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER} \
    ${GAPS_COORDINATES_LG} \
    ${SCAFFOLD_FA}
# Output: $SCAFFOLD_FA - scaffold sequences with Ns in gaps

echo "---> Computing ancestral sequence with gaps"
python2 src/fpsac_compute_scaffold_sequence.py \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER} \
    ${GAPS_COORDINATES_LG} \
    ${FINISHING_ALG}/ \
    ${SCAFFOLD_FINISHED_FA}
# Output: $SCAFFOLD_FINISHED_FA - FASTA sequence of the ancestral scaffold, including gaps

echo "---> Computing ancestral map"
python2 src/fpsac_compute_scaffold_map.py \
    ${FAMILIES_COORDINATES_ALL} \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER} \
    ${GAPS_COORDINATES_LG} \
    ${FINISHING_ALG}/ \
    ${SCAFFOLD_MAP}
# Output: $SCAFFOLD_MAP - map of the ancestral scaffold

# -----------------------------------------------------------
# echo
# echo "STAGE D: GENERATING REPORT"
# 
# python src/fpsac_generate_report.py \
#     $REPORT \
#     $ANCESTRAL_CONTIGS_LG \
#     $EXTANT_GENOMES_LG \
#     $ANCIENT_EXTANT_HITS \
#     $FAMILIES_COORDINATES_ALL \
#     $FAMILIES_ANCESTRAL_CONTENT \
#     $FAMILIES_COORDINATES_FILTERED \
#     $FAMILIES_LG \
#     $ADJACENCIES_WGTD \
#     $ADJACENCIES_KEPT \
#     $ADJACENCIES_DISC \
#     $REPEAT_CLUSTERS \
#     $RSI_WGTD \
#     $RSI_KEPT \
#     $RSI_DISC \
#     $ADJACENCIES_CARS_ALL \
#     $SCAFFOLD_ORDER \
#     $GAPS_COORDINATES_LG \
#     $SCAFFOLD_MAP \
#     $ADJACENCIES_UNASSIGNED


