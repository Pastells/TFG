#!/bin/bash
echo $(date)
main() {
    local rate=$1
    for lambda in 20
        do
        echo $rate, $lambda
        mkdir temp_${rate}_${lambda}
        cp sir_mut temp_${rate}_${lambda}
        cp array.dat temp_${rate}_${lambda}
        cp ../inp/sir_mut.inp temp_${rate}_${lambda}
        cd temp_${rate}_${lambda}
        sed -i "s/mutation_rate=.*/mutation_rate=${rate}d${Order},/" sir_mut.inp
        sed -i "s/lambda=.*/lambda=${lambda}d0,/" sir_mut.inp
        ./sir_mut
        cd ..
        cat temp_${rate}_${lambda}/endemic.dat >> endemic.dat
    done
}

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
FOLDER="sir_mut_$DATE"
mkdir $FOLDER
cp edge_list.dat $FOLDER
cp degree.dat $FOLDER

echo 'rate, lambda'
Order=2
for rate in {732..1326..66}
    do main "$rate" &
done
wait
mv endemic.dat $FOLDER
rm -rf temp*
cd ~/tfg/simulations
mv programs/bin/$FOLDER results/$FOLDER
echo $(date)
