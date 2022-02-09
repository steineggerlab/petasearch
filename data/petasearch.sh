#!/bin/sh -e
# petasearch workflow

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}


#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 5 ] && echo "Please provide <queryDB> <targetDBsFile> <compResult> <alignmentRes> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
#[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -d "$5" ] && echo "tmp directory $5 not found! Creating it." && mkdir -p "$5";


Q_DB="$1"
T_DBs="$2"
C_RES="$3"
SA_RES="swappedAlis.out"
ALI_RES="alis.out"
M8_RES="alis.m8"
FINAL_RES="$4"
TMP_PATH="$5"




# compare both k-mer tables
#if [ ! -e "${TMP_PATH}/${C_RES}"  ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" compare2kmertables "${Q_DB}" "${T_DBs}" "${C_RES}" ${COMP_KMER_TABLES_PAR} \
        || fail "comparing both k-mer tables failed"
#fi
paste "${T_DBs}" "${C_RES}"  | column -s "\t" > "threecol.tsv"

STEP=0
# shellcheck disable=SC2034
while IFS="$(printf "\t")" read -r TARGETABLE T_DB COMP_RES; do
   # shellcheck disable=SC2086
  "$MMSEQS" computeAlignments "${Q_DB}" "${T_DB}" "${COMP_RES}" "${TMP_PATH}/${ALI_RES}_${STEP}" ${COMP_ALI_PAR} \
        || fail "computing the alignment for matched sequences failed"

   # shellcheck disable=SC2086
#  "$MMSEQS" swapresults "${T_DB}" "${Q_DB}" "${TMP_PATH}/${SA_RES}_${STEP}" "${TMP_PATH}/${ALI_RES}_${STEP}" ${SWAP_PAR} \
#        || fail "swapping the alignments failed"
#
#  # shellcheck disable=SC2086
#  "$MMSEQS" convertalis "${Q_DB}" "${T_DB}" "${TMP_PATH}/${ALI_RES}_${STEP}" "${TMP_PATH}/${M8_RES}_${STEP}" ${CONVERTALIS_PAR} \
#      || fail "creating  the .m8 file failed"

  # shellcheck disable=SC2086
  "$MMSEQS" convertsraalis "${Q_DB}" "${T_DB}" "${TMP_PATH}/${ALI_RES}_${STEP}" "${TMP_PATH}/${M8_RES}_${STEP}" 
  ${CONVERTALIS_PAR} \
      || fail "creating  the .m8 file failed"

  push_back "${TMP_PATH}/${M8_RES}_${STEP}"
  STEP=$((STEP+1))
done < "threecol.tsv"

eval "set -- $ARR"
cat "${@}" > "${FINAL_RES}"

# clear up tmp files
if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}"
fi



