#!/bin/sh
	
echo "Please enter your num1"
echo "followed by num 2 path and name eg 30 49 ./transposed_tables/ multiplex_  \c"
read START END PAT MP

   
for i in `seq $START $END`
do
	#echo $i
	cmd="mkdir ${PAT}${MP}${i}"
	echo $cmd
	$cmd
	
done

# ./scripts/make_dir.sh
# 1 26 ./recoded_loc_files/ multiplex_
# 1 26 ./transposed_tables/ multiplex_
# 1 26 ./single_loc_files/ multiplex_
