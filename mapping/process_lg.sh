#!/bin/sh
echo "Please enter the file that contains the LG names"
echo "followed by any extra loci to append   \c"
cd $PWD
read p LOC EXT

#loc_from_list.pl lg3b.txt all.loc lg3b_list.loc extra.txt


while read p; do
	cmd="loc_from_list.pl ${p}.txt $LOC ${p}_list.loc ${EXT} ${p}"
  	echo $cmd
  	$cmd
done < $p
