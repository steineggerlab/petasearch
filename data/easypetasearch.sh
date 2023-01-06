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

# $1: query database; $2: target database; $3 alignment result; $4 tmp
TMP_PATH="$4"

if [ ! -e "${TMP_PATH}/res" ]; then
  mkdir -p "${TMP_PATH}/res"
fi

if [ ! -e "${TMP_PATH}/tmp" ]; then
  mkdir -p "${TMP_PATH}/tmp"
fi

  # shellcheck disable=SC2086
  "$MMSEQS" petasearch "$1" "$2" "${TMP_PATH}/res" "$3" "${TMP_PATH}/tmp" ${PETASEARCH_PAR} \
        || fail "petasearch failed"

# This procedure has been added to petasearch workflow
# shellcheck disable=SC2086
#"$MMSEQS" convertalis "$1" "$2" "${TMP_PATH}/res" "$3" ${CONVERTALIS_PAR} \
#      || fail "creating  the .m8 file failed"

# clear up tmp files
if [ -n "$REMOVE_TMP" ]; then
    echo "Remove temporary files"
    "$MMSEQS" rmdb "${TMP_PATH}/res"
    rm -rf "${TMP_PATH}"
fi
