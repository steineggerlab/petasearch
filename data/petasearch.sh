#!/bin/sh -e
# petasearch workflow

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

post_proc () {
    STEP="$1"
    # shellcheck disable=SC2086
    "$MMSEQS" blockalign "${Q_DB}" "${T_DB}" "${COMP_RES}" "${TMP_PATH}/${ALI_RES}_${STEP}" ${COMP_ALI_PAR} \
        || fail "computing the alignment for matched sequences failed"

    # shellcheck disable=SC2086
    "$MMSEQS" convertsraalis "${Q_DB}" "${T_DB}" "${TMP_PATH}/${ALI_RES}_${STEP}" "${TMP_PATH}/${M8_RES}_${STEP}" ${CONVERTALIS_PAR} \
        || fail "creating  the .m8 file failed"
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 5 ] && echo "Please provide <queryDB> <targetDBsFile> <compResult> <alignmentRes> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
#[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;

Q_DB="$1"
T_DBs="$2"
C_RES="$3"
#SA_RES="swappedAlis.out" # Suppress SC2034
ALI_RES="alis.out"
M8_RES="alis.m8"
FINAL_RES="$4"
TMP_PATH="$5"

# compare both k-mer tables
if notExists "${TMP_PATH}/comparekmertables.done"; then
    # shellcheck disable=SC2086
    "$MMSEQS" comparekmertables "${Q_DB}" "${T_DBs}" "${C_RES}" ${COMP_KMER_TABLES_PAR} \
        || fail "comparing k-mer tables failed"
    touch "${TMP_PATH}/comparekmertables.done"
fi

paste "${T_DBs}" "${C_RES}"  | column -s "\t" > "${TMP_PATH}/threecol.tsv"

STEP=0
# shellcheck disable=SC2034
while IFS="$(printf "\t")" read -r TARGETABLE T_DB COMP_RES; do
    # shellcheck disable=SC2086
    post_proc $STEP &
#  push_back "${TMP_PATH}/${M8_RES}_${STEP}"
    STEP=$((STEP+1))
done < "${TMP_PATH}/threecol.tsv"
wait

STEP=$((STEP-1))
: > "${FINAL_RES}" # Suppress SC2188
for i in $(seq 0 $STEP); do
    cat "${TMP_PATH}/${M8_RES}_${i}" >> "${FINAL_RES}"
done

# # clear up tmp files
# if [ -n "$REMOVE_TMP" ]; then
#     echo "Remove temporary files"
#     rm -f "${TMP_PATH}/threecol.tsv"
#     rm -rf "${TMP_PATH}"
# fi
