#! /bin/bash -e

if (( $# < 3 )); then
    cat >&2 <<EOM
usage: $0 <introns+> <introns-> <fasta>
    introns+: positive introns dump
    introns-: negative introns dump
    fasta: fasta file (having generated an index for it)
EOM
    exit 1
fi

echo "Parsing positives..."
while read LINE; do
    #echo "${LINE##*,}..."
    >>parsedpos.txt samtools faidx "$3" "${LINE##*,}"
done < "$1"

echo "Parsing negatives..."
while read LINE; do
    #echo "${LINE##*,}..."
    >>parsedneg.txt samtools faidx -i "$3" "${LINE##*,}"
done < "$2"
