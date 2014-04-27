#!/bin/sh



cd  /$HOME/projects/gene_identification/results/blast_results/

for i in `ls bla*`
do 
cd /$HOME/projects/gene_identification/
echo $i 
./term_finder.pl ./results/blast_results/$i ./data/ATH_GO_GOSLIM.txt  ./data/keywords.txt ./data/go_terms.txt ./results/filtered/$i
done

