#!/bin/bash

# Should be executed as: ./all_sis.sh input.dat
NAME=${1?Error: no filename given}
#cat $NAME | sort -n | uniq -u > programs/bin/edge_list.dat
sort -n $NAME | uniq -u > programs/bin/edge_list.dat
echo 'edge_list.dat created'

cd programs/bin
# Check if enumeration starts with zero or one
read -r FIRSTLINE < edge_list.dat
if [[ $FIRSTLINE == 0* ]]; then ./read_zero; else ./read_one; fi
echo 'array.dat created'

DATE=$(date +%F_%H_%M)
FOLDER="sis_$DATE"
mkdir $FOLDER
cp edge_list.dat $FOLDER
cp degree.dat $FOLDER

./sis > $FOLDER/prints.dat
mv mean.dat $FOLDER/mean.dat
cd -
mv programs/bin/$FOLDER results/$FOLDER
