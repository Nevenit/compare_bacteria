iters=100
start=`date +%s.%N`
for i in `seq 1 $iters`
do
./target/release/compare_bacteria list.txt
done
end=`date +%s.%N`

runtime=$( echo "($end - $start) / $iters" | bc -l )
echo "The average time of $iters runs is: $runtime s"
