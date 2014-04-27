
#blast the vesca genes against custom blast db's of the scaffold and retrieve the genes on those scaffolds...
#IN THE SCRIPTS DIRECTORY

./vescablaster_gff.pl /home/harrir/projects/gene_identification/data/scaffolds.txt /home/harrir/projects/gene_identification/results/map_gff.txt

#find genes homologous to vesca genes on the scaffold
#IN THE GENE IDENTIFICATION DIRECTORY

./vesca_annotator.pl ./results/map_gff.txt ./data/strawberry_gene_model_plus_dna_161109.fna TAIR10_cdna_20101214_NEW.txt ./results/blast_results/blast


#filter the results on keywords present in homologs and on ontology terms that arab homologs contain. This puts the data in the ./results/filtered directory
#IN THE SCRIPTS DIRECTORY

./go_filter.sh 

#this runs this script (but for all files)
#./term_finder.pl ./results/blast_results/reblast_scf0512956.txt ./data/ATH_GO_GOSLIM.txt  ./data/keywords.txt ./data/go_terms.txt ./results/filtered/scf0512956_filtered.txt

#MANUALLY select candidates from blast and ontology and place scaffold file with list of mrna's into the candidate directory

# design primers of selected candidates in introns where possible
#IN THE SCRIPTS DIRECTORY

./bulk_primer_v2.pl /home/harrir/projects/gene_identification/results/candidates/ /home/harrir/projects/gene_identification/results/primers/ /home/harrir/projects/gene_identification/results/seq/ /home/harrir/projects/gene_identification/results/summary/


#get some info about the specific genes chosen by the candidate search
#IN THE SCRIPTS DIRECTORY

#./candidate_blast.pl /home/harrir/projects/gene_identification/results/candidates/ /home/harrir/projects/gene_identification/results/map_candidates.txt

#IN THE GENE IDENTIFICATION DIRECTORY

#./vesca_annotator.pl ./results/map_candidates.txt ./data/strawberry_gene_model_plus_dna_161109.fna TAIR10_cdna_20101214_NEW.txt ./results/reblast_results/reblast

