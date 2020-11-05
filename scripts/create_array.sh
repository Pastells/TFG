#!/bin/bash

# Should be executed as: ./create_array.sh input.dat
NAME=${1?Error: no filename given}
#cat $NAME | sort -n | uniq -u > programs/bin/edge_list.dat
sort -n $NAME | uniq -u > programs/bin/edge_list.dat
echo 'edge_list.dat created'

cd programs/bin
# Check if enumeration starts with zero or one
read -r FIRSTLINE < edge_list.dat
if [[ $FIRSTLINE == 0* ]]; then ./read_zero; else ./read_one; fi
echo 'array.dat created'
echo 'degree.dat created'
