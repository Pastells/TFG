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


for lambda in {55..100..5}
    do sed -i "s/lambda=.*/lambda=${lambda}d0,/" ../inp/sir_mut_i.inp
    echo ${lambda}
    ./sir_mut_i
    #./sir_mut > $FOLDER/rate${rate}.out
done
cd -
