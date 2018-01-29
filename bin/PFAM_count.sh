#!/bin/bash
# Count number of sequences of each fasta file and move them to a folder PDB_10seqs

echo -n "Enter the minimun number of sequences for each PFAM family > "
read size
echo -n "Create a folder to store the sequence files > "
read folder
    mkdir $folder
echo -n "Create a folder to store unused sequence files > "
read folder2
    mkdir $folder2

for entry in `ls $PWD/*pdb.fa`; do
    seqs=$(grep \> $entry | wc -l) 
    if [ "$seqs" -ge $size ]; then
    fam=$(echo "Family $entry" | awk -F/ '{print $NF}' | cut -f1 -d "_")
    files="${fam}*"
    echo $files
    mv $files $folder
    fi
    if [ "$seqs" -lt $size ]; then
    fam=$(echo "Family $entry" | awk -F/ '{print $NF}' | cut -f1 -d "_")
    files="${fam}*"
    echo $files
    mv $files $folder2
    fi

done
echo $size
