#!/bin/bash

# Should be executed as: ./all_sis_endemic.sh input.dat
NAME=${1?Error: no filename given}
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
touch endemic_sis.dat

# lambdas to loop over
    #L0=150
    #Lf=220
    #Lstep=10
#for (( lambda=$L0; lambda<=$Lf; lambda+=$Lstep ))
for lambda in {150..180..1}
    do sed -i "s/lambda=.*/lambda=${lambda}d0,/" ../inp/sis_endemic.inp
    echo $lambda
    ./sis_endemic
done
cp endemic_sis.dat $FOLDER
echo >> endemic_sis.dat
echo >> endemic_sis.dat
cd -
mv programs/bin/$FOLDER results/$FOLDER
