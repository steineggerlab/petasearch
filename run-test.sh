#!/usr/bin/zsh
# run this only after build

BRANCH=''
k=''

# b: switch between branches; k: kmer length
while getopts b:k: option
do
    case "${option}" in
      b) BRANCH=${OPTARG};;
      k) k=${OPTARG};;
      *) ;;
    esac
done

currbr=$(git rev-parse --abbrev-ref HEAD)

if [ -n "$BRANCH" ]; then
    git checkout "$BRANCH"
fi

br=$(git rev-parse --abbrev-ref HEAD)

if [ -z "$k" ]; then
    k=11
fi

cd ~/summer-research/SRA-search/build/ || exit
make -j 16
cd ..

# Treat vanilla (master branch result) specially
if [[ $br == "master" ]] ; then
    res="vanilla(k=${k}).log"
else
    res="analysis.${br}(k=${k}).log"
fi

# It seems that we don't have to do this
start=$(date +%s)
~/summer-research/SRA-search/build/src/srasearch createkmertable /home/matchy233/summer-research/database/swissprot /home/matchy233/summer-research/database/swissprot_kmertable -k $k
end=$(date +%s)
runtime=$((end-start))

hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))

echo "### createkmertable runtime ###" > "$res"
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)" >> "$res"
echo " " >> "$res"

start=$(date +%s)
~/summer-research/SRA-search/build/src/srasearch compare2kmertables --exact-kmer-matching 1 /home/matchy233/summer-research/database/swissprot /home/matchy233/summer-research/target_table_list /home/matchy233/summer-research/result_output_path -k $k
end=$(date +%s)
runtime=$((end-start))

hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))

echo "### compare2kmertables runtime ###" >> "$res"
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)" >> "$res"
echo " " >> "$res"

if [[ $br == "master" ]] ; then
    out="vanilla(k=${k}).tsv"
    if [ ! -f "~/summer-research/vanilla(k=${k}).tsv" ] ; then
        echo > "$out"
    fi
    ~/summer-research/SRA-search/build/src/srasearch createtsv ~/summer-research/database/swissprot ~/summer-research/database/swissprot ~/summer-research/result/result_swissprot $out
    mv $out ~/summer-research/
else
    ~/summer-research/SRA-search/build/src/srasearch createtsv ~/summer-research/database/swissprot ~/summer-research/database/swissprot ~/summer-research/result/result_swissprot ~/summer-research/tmp.tsv
    FILE=/home/matchy233/summer-research/vanilla\(k=${k}\).tsv
    if [ ! -f "$FILE" ] ; then
        echo "No vanilla result to check the correctness!"
        exit 1
    fi

    FILE=/home/matchy233/summer-research/vanilla\(k=${k}\).tsv
    curr=$(wc -c < ~/summer-research/tmp.tsv)
    old=$(wc -c < "$FILE")

    if [ "$curr" -ne "$old" ]; then
        echo "Current implementation is wrong!"
        exit 1
    fi
fi

kmerdb=$(wc -c < ~/summer-research/database/swissprot_kmertable)
iddb=$(wc -c < ~/summer-research/database/swissprot_kmertable_ids)

echo "=== kmer table size ===" >> "$res"
echo "$kmerdb" >> "$res"
echo " " >> "$res"

echo "===kmer id size ===" >> "$res"
echo "$iddb" >> "$res"

rm -f ~/summer-research/tmp.tsv
mv "$res" ~/summer-research/

if [ -n "$BRANCH" ]; then
    git checkout "$currbr"
fi
