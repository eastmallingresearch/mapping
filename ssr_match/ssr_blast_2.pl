#!/usr/bin/perl


use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqFeature::Generic;
use Bio::LocatableSeq;
use Bio::DB::SeqFeature;
use Bio::DB::SeqFeature::Store;


# path to text file of contigs
 my $contig = shift;
# path to output 
 my $output = shift;

#./ssr_blast.pl ssr.csv >output.txt
#formatdb -p F -i fvesca.090806.scf.fna -n vescadb


#array of contigs
print "Collecting contigs \n";
open(IN, $contig); #open for write
my @contigs=<IN>;
close IN;


my %contig_gene=();
# print "Genes being being examined\n";
	
	foreach (@contigs){
		chomp;
		chop;
	#	print $_."\n";
		my @keyval=split(",",$_);
		#print $keyval[0]."\n";
	#my $seq=();
		my $seq = Bio::Seq->new( -seq => $keyval[1],
                                 -id  => $keyval[0],
				 
			       );
			       
		$contig_gene{$keyval[0]}=$seq;	
		
	}
#exit;
	
foreach (keys %contig_gene){
		#print "Working on primer ".$_."\n";
		#print "Conducting comparison to vesca genome\n";
		#my %hash=%{$simplified_hash{$_}};
		my @hits=blaster(\$contig_gene{$_});
		#my $filename=$output."_".$_.".txt";
		#print "Outputting Data for scaffold $_ \n";
		#blast_output(\%hits,$filename,$_);
		#$hash_of_results{$_}=\%hits;		
		}
	





sub blaster{
	
	my ($seq)=@_;
	
	my @params = (program  => 'blastn', database => "./blast_scaff/vescadb" );
	my	$blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
	 my @aohits=();
       			#print $$seq->id()."\n";
       	 		my $query_length= $$seq->length();
   			#	print "Query Length $query_length\n";
       			my  $report_obj = $blast_obj->blastall($$seq);
       			my %hits=report_filter($report_obj,$$seq);        					
       			#print keys %hits;
       			
       			if (scalar(keys %hits)>0){
       			##	print ("RETURNING VALID HITS\n");
       				push (@aohits,\%hits);

				foreach( keys %hits){
					print ($hits{$_}->id(),"\t",$_,"\n");
						}
				#exit;
       			}
       	#exit;
       	return (@aohits);
	
}

sub report_filter{
	my ($report_obj,$seq)=@_;
	my $query_length=$seq->length();
	my %hashofhits=();
	#print ("Query Length ", $query_length,"\n");
		while(my $result = $report_obj->next_result ) {
       				my $hitcount= $result->num_hits;
				print "HITS $hitcount \n";
 			 			 while( my $hit = $result->next_hit ) {
								print $hit->name()."\t";
   							while( my $hsp = $hit->next_hsp ) {
							    print $hsp->query_string."\t";
									print $hsp->hit_string."\t";
							                print $hsp->length('query')."\t";
									print $hsp->start('hit')."\t";
									print $hsp->end('hit')."\t";
							    print $hsp->strand('hit')."\t";
							                print $hsp->percent_identity."\n";


          								if( $hsp->length('total') == $query_length && $hsp->frac_identical( ['query'|'hit'|'total']) >0.9  ) {
#if( $hsp->length('total') == $hsp->frac_identical( ['query'|'hit'|'total']) >0.5  ) {

     												$hashofhits{$hit->name()}=$seq;
											print $hit->name()."\n";
     											print "Significant hit on scaffold \n";
     											        	
          								}
    				       			}
			 			 		}
       			 			}

		
		return %hashofhits;	
}

