#!/usr/bin/perl


use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;


# path to text file of results
 my $results = shift;
# path to gene ontology
 my $go = shift;
# path to gene keywords
 my $keypath = shift;
 # path to go keywords
 my $gokey = shift;
# path to output 
 my $output = shift;

#./term_finder.pl ./results/reblast_scf0512956.txt ./data/ATH_GO_GOSLIM.txt  ./data/keywords.txt ./data/go_terms.txt scf0512956_filtered.txt


#read in data and keywords
my @raw_results=read_in($results);
my @godata=read_in($go);
my @keywords=read_in($keypath);
my @go_keywords=read_in($gokey);

#print @raw_results;
#print @godata;
#print @keywords;
#print @go_keywords;
#exit;

#process the raw data 
my %gbg=genebygene(@raw_results);

#print keys %gbg;
#exit;

print "Keyword Scores\n";
#get keyword scores
my %keyword_scores=keyword_scores(\%gbg,\@keywords);
print "Ontology Scores\n";
#process the go_data into a hash
my%go_proc=process_ontology(\@godata);
my %ontology_scores=ontology_scores(\%gbg, \@go_keywords,\%go_proc);

output_scores(\%keyword_scores,\%ontology_scores,$output);


#foreach (keys %keyword_scores){
#		print $_."\t";
#		my @tmp=@{$keyword_scores{$_}};
#	
#		print $tmp[0],"\t",$tmp[1],"\t";
#			foreach (@{$tmp[2]}){
#				print $_."\t"
#			}
#			print "\n";
#	
#		}
#
#
#foreach (keys %ontology_scores){
#		print $_."\t";
#		my @tmp=@{$ontology_scores{$_}};
#	
#		print $tmp[0],"\t",$tmp[1],"\t";
#			foreach (@{$tmp[2]}){
#				print $_."\t"
#			}
#			print "\n";
#	
#		}

##########################SUBROUTINES#############################

sub output_scores{
my ($key,$ont,$out)=@_;

my %keyword=%{$key};
my %ontology=%{$ont};
print "OUTPUTTING TO $out \n";

open(OUT, ">$out"); #open for read



foreach my $key (keys %keyword){
	foreach (keys %ontology){
		if ($key eq $_){
		
			my @keys=@{$keyword{$key}};
			my @onts=@{$ontology{$key}};	
			
			
			
			print $key,"\t",$keys[0],"\t",$keys[1],"\t";
			print OUT $key,"\t",$keys[0],"\t",$keys[1],"\t";
			
			foreach (@{$keys[2]}){
				print $_.",";
				print OUT $_.",";
			}
			
			print "\t",$onts[0],"\t",$onts[1],"\t";
			print OUT "\t",$onts[0],"\t",$onts[1],"\t";
			
			foreach (@{$onts[2]}){
				print $_.",";
				print OUT $_.",";
			}
			print "\n";
			print OUT"\n";

		}#if
		}#foreach
	}#foreach
close (OUT);

}


sub read_in{
my ($path)=@_;
print "Collecting results into groups $path\n";

open(IN, $path); #open for read
my @raw_results=<IN>;
close IN;
return @raw_results;
}



sub process_ontology{
my ($go)=@_;
my @raw_results=@$go;
my %ho_genearrays=();
my @array;

	for (my $i=1;$i<@raw_results;$i++){
			my @temp_array_b=split("\t",$raw_results[$i-1]);
			my @temp_array_f=split("\t",$raw_results[$i]);
			
				#print $temp_array_b[1]."\n";
				
				if ($temp_array_b[0] eq $temp_array_f[0]){
					
					push(@array,$temp_array_b[4]);
				}
				else{
					#print "adding ",$temp_array_f[1]."\n";
					my @new_array=@array;	
					
					#	print $temp_array_b[0]."\n";
					#	foreach (@new_array){
					#		print $_."\n";
					#	}
						#exit;
						
					$ho_genearrays{$temp_array_b[0]}=\@new_array;
					@array=();
				} 
		}
	return %ho_genearrays;
	
	


}


sub keyword_scores{
my ($gb,$key)=@_;
my %keyscores=();
my %gbg=%{$gb};

	foreach(keys %gbg){
	
		print "DATA FOR $_ \n";
		my $score=0;
		my %keywords=();
		my @array=@{$gbg{$_}};
		my $length=scalar(@array);
		
		foreach (@array){
			my @orf=split("\t",$_);
			#$orf[2] =~ tr/[A-Z]/[a-z]/;
			#$orf[3] =~ tr/[A-Z]/[a-z]/;
			#$orf[4] =~ tr/[A-Z]/[a-z]/;
			#print $orf[2],$orf[3],$orf[4],"\n";
				foreach (@$key){
					#print $_, "\t",$orf[3],"\n";
					chomp $_;
					if( $orf[2]=~m/$_/i || $orf[3]=~m/$_/i || $orf[4]=~m/$_/i){
						#print "MATCHING TERM $_ \n";
						$score=$score+1;
						$keywords{$_}++;
					}
				}
		}
		my @score_words=();
		push (@score_words,$score);
		push (@score_words,$length);
		my @tmp=keys(%keywords);
		push (@score_words,\@tmp);
		$keyscores{$_}=\@score_words;
		#print "Keyword score $_ $score/$length ", keys (%keywords)," \n";
	}
return %keyscores;
}



sub match_thal{

	my ($gene,$thal)=@_;
	my @aot=();
	my @matches=();
	my %data=%{$thal};
	
	#collect thal refs for a gene
	foreach ($gene){
		my @tmp=@{$_};
		foreach (@tmp){
			my @splits=split ('\t', $_);
			my @sub_split=split('_',$splits[2]);
			 push(@aot,$sub_split[0]);
		}
	}

	#compare thal refs to find matches- return an array
	foreach my $goref (keys %data){
		 foreach (@aot){
	 		if ($_=~/$goref/){
	 		#print "MATCH $_ $goref \n";
	 	
	 		push (@matches,$data{$goref});
	 		}
		 }
	}

	my @revealed=();

	foreach (@matches){
		foreach (@{$_}){
			push(@revealed, $_);
		}
	}

return (@revealed);


}


sub ontology_scores{

my ($gb,$key,$god)=@_;
my %goscores=();
my %gbg=%{$gb};
my %godata=%{$god};

foreach(keys %gbg){
	

		
	
		#print "DATA FOR $_ \n";
		my $score=0;
		my %keywords=();
		my @array=@{$gbg{$_}};
		my @matches=match_thal(\@array,\%godata);
		
		my $length=scalar(@array);
		foreach my $match (@matches){
			foreach (@$key){
				#print $match,"\t",$_,"\n";
					chomp $_;
					if( $match=~m/$_/  ){
						#print "MATCHING TERM $_ \n";
						$score=$score+1;
						$keywords{$_}++;
					}
				}
		}
		my @score_words=();
		push (@score_words,$score);
		push (@score_words,$length);
		my @tmp=keys(%keywords);
		push (@score_words,\@tmp);
		$goscores{$_}=\@score_words;
		#print "Keyword score $_ $score/$length ", keys (%keywords)," \n";
		}
	return %goscores;
	
}

sub genebygene{

	my ($data)=@_;
	my %ho_genearrays=();
	my @array=();			
				
		for (my $i=2;$i<@raw_results;$i++){
			my @temp_array_b=split("\t",$raw_results[$i-1]);
			my @temp_array_f=split("\t",$raw_results[$i]);
							
				if ($temp_array_b[1] eq $temp_array_f[1] && $i !=(@raw_results-1) ){
				#print "IN MAIN";
					push(@array,$raw_results[$i-1]);
					
				}
								
				else{
				#print "HERE";
					#print "adding ",$temp_array_f[1]."\n";
					my @new_array=@array;	
					$ho_genearrays{$temp_array_b[1]}=\@new_array;
					@array=();
				} 
		}
		
		
	return %ho_genearrays;
}



