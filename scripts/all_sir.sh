#!/bin/bash

# lambdas to loop over
    L0=10
    Lf=200
    Lstep=5

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

# total infected
    sed -i "s/lambda_0=.*/lambda_0=${L0},/" total_infected_sir.inp
    sed -i "s/lambda_f=.*/lambda_f=${Lf},/" total_infected_sir.inp
    sed -i "s/lambda_step=.*/lambda_step=${Lstep},/" total_infected_sir.inp
    ./total_infected_sir


DATE=$(date +%F_%H_%M)
FOLDER="sir_$DATE"
mkdir $FOLDER
cp edge_list.dat $FOLDER
cp degree.dat $FOLDER
cp total_inf.dat $FOLDER
touch recovered.dat
touch survival.dat

for (( lambda=$L0; lambda<=$Lf; lambda+=$Lstep ))
#for lambda in {${L0}...${Lf}...${Lstep}}
    do sed -i "s/lambda=.*/lambda=${lambda}d0,/" ../inp/sir.inp
    #./sir > $FOLDER/lambda_${lambda}.out
    echo $lambda
    ./sir
done
mv recovered.dat $FOLDER/recovered.dat
mv survival.dat $FOLDER/survival.dat
cd -
mv programs/bin/$FOLDER results/$FOLDER
