#!/bin/sh -e
# petasearch workflow

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <alignmentRes> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found! Creating it." && mkdir -p "$4";


Q_DB="$1"
T_DB="$2"
T_TABLE="targetTable"
C_RES="compResults.out"
SA_RES="swappedAlis.out"
FINAL_RES="$3"
TMP_PATH="$4"


#create target table
if [ ! -e "${TMP_PATH}/${T_TABLE}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" createkmertable "${T_DB}" "${TMP_PATH}/${T_TABLE}" ${CREATE_TTABLE_PAR} \
        || fail "creating target table failed"
fi

# compare both k-mer tables
if [ ! -e "${TMP_PATH}/${C_RES}"  ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" compare2kmertables "${Q_DB}" "${TMP_PATH}/${T_TABLE}" "${TMP_PATH}/${C_RES}" ${COMP_KMER_TABLES_PAR} \
        || fail "comparing both k-mer tables failed"
fi

# compute the alignment for significant matches
if [ ! -e "${TMP_PATH}/${SA_RES}"  ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" computeAlignments "${Q_DB}" "${T_DB}" "${TMP_PATH}/${C_RES}" "${TMP_PATH}/${SA_RES}" ${COMP_ALI_PAR} \
        || fail "computing the alignment for matched sequences failed"

fi

# swap the alignments
if [ ! -e "${FINAL_RES}"  ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapresults "${T_DB}" "${Q_DB}" "${TMP_PATH}/${SA_RES}" "${FINAL_RES}" ${SWAP_PAR} \
        || fail "swapping the alignments failed"
fi

# clear up tmp files
if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}"
fi



