mapping
=======

Linkage mapping and primer design

#place a list of primers into a text file named multiplex_n in the directory ./data_analysis/multiplex_names/raw. This script goes through the data 
and removes those that have failed from the multiplexes and those that have not yet been scored from the multiplexes and appends them to the 
unscorable list. #ret 0 = normal 1=fail  2=unscore 

#./scripts/generate_lists.pl ./edited_tables/checked_csv/ ./edited_tables/unscorable_list/ ./multiplex_names/raw/ ./multiplex_names/processed/ 27 49

#./scripts/generate_lists.pl ./edited_tables/checked_csv/ ./edited_tables/unscorable_list/ ./multiplex_names/raw/ ./multiplex_names/processed/ 9 9


./scripts/make_dir.sh 1 49 ./recoded_loc_files/ multiplex_


#create a directory named multiplex_n in the ./transposed_tables directory and the ./single_loc_files directory. This script will then call the marker
type as nnxnp, lmxll or hkxhk and output the marekr into the transposed_tables and single_loc files directory. This then needs manually scoring in joinmap

#./scripts/format_col.pl ./edited_tables/checked_csv/ ./multiplex_names/processed/ ./transposed_tables/ ./single_loc_files/ multiplex_27
#./scripts/format_col.pl ./edited_tables/checked_csv/ ./multiplex_names/processed/ ./transposed_tables/ ./single_loc_files/ test



