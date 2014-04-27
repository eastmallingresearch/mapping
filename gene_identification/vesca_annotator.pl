#!/usr/bin/perl


use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;

# path to text file of contigs
 my $contig = shift;
# path to genome sequence
 my $genome_mrna = shift;
# path to TAIR file 
 my $compgen = shift;
# path to output 
 my $output = shift;

#./vesca_annotator.pl ./results/newrh_genes.txt ./data/gene_model.txt TAIR10_cdna_20101214_NEW.txt reblast



#array of contigs
print "Collecting contigs \n";
open(IN, $contig); #open for write
my @contigs=<IN>;
close IN;



#exit;
my %contig_gene=();

print "Genes being being examined\n";
	
	foreach (@contigs){
		chomp;
		#print $_."\n";
		my @keyval=split("\t",$_);
		#print $keyval[0]."\n";
		$contig_gene{$keyval[1]}=$keyval[0];	
		
	}
	
#print keys %contig_gene;



#exit;
print "Searching Genome $genome_mrna \n";

my %genes=search_genome(\%contig_gene,$genome_mrna);

#print "

print "Scaffolds identified and sequences stored in text file and blast database \n";
print "Simplifying hash\n";

my %simplified_hash=simplify_hash(\%genes,\@contigs);


#comparatives 

	foreach (keys %simplified_hash){
		print "Working on scaffold ".$_."\n";
		print "Conducting comparison to model organism\n";
		my %hash=%{$simplified_hash{$_}};
		my %hits=arab_blaster(\%hash,\$compgen);
		my $filename=$output."_".$_.".txt";
		print "Outputting Data for scaffold $_ \n";
		blast_output(\%hits,$filename,$_);
		#$hash_of_results{$_}=\%hits;		
		}
		
	
	
	
	print "BLAST complete\n";

#####SUBROUTINES 

sub simplify_hash{
	my($genes,$list)=@_;
	
	my %hashofscaffolds;
	my %genes=%{$genes};
	my @list=@$list;
	print "matching gene to scaffolds\n";
	my %scaffolds;
	
	foreach (@list){
		chomp;
		print $_."\n";
		my @keyval=split("\t",$_);
		$scaffolds{$keyval[0]}=$keyval[1];	
	}
	
	
	foreach my $scaffold (keys %scaffolds){
		my %scaffold_hash=();
		print "Working on scaffold $scaffold \n";
			
			for (my $i=0; $i<scalar(@list);$i++){
				my @keyval=split("\t",$list[$i]);
				print "Working on gene_id ".$keyval[1]." and searching for it \n";
			
					if ($keyval[0] eq $scaffold){
						foreach my $orf (keys %genes){
							my $seq=$genes{$orf};
								if ( $seq->id() eq $keyval[1]){
									print "We have a match $keyval[0] ". $seq->id()."\t". $orf." \n";
									$scaffold_hash{$keyval[1]}=$seq;
								}	
							}
					}
					else{
						next;
					}
				}
		
		print "adding hash of genes in scaffold to hash of scaffolds\n";
		$hashofscaffolds{$scaffold}=\%scaffold_hash;
	}
	
	
	return %hashofscaffolds;
}


sub search_genome{
	my($contigs,$genomepath )=@_;
	
	
	my %hash=();
#	print keys %{$contigs};
	
	my $inseq = Bio::SeqIO->new(-file => "<$genomepath",
                                -format => 'fasta' );
                                
    
      
       while (my $seq = $inseq->next_seq) {
       		
         	foreach(keys %{$contigs}){
        #		chomp $_;
         		
         			if ($seq->id eq $_){
         				print "Found ". $seq->id()."\n"; 
         				$hash{$_}=$seq;		
         			}
         		}
    		}


	return %hash;
	
}

sub create_database{
	my ($contigs)=@_;
	print "Creating Local Blast Database \n";
	my $seq_out = Bio::SeqIO->new('-file' => ">out_tmp.fa",
                                       '-format' => 'fasta');
	my %hash=%{$contigs};
	
	foreach (keys %hash ){
	
	my $seq_out = Bio::SeqIO->new('-file' => ">>out_tmp.fa",
                                       '-format' => 'fasta');

         # write each entry in the input file to the output file
        	
            $seq_out->write_seq($hash{$_});
	}
         
	#print $contig->seq();
	my @argv=("./temp_blast.sh","out_tmp.fa");
	system(@argv);




}

sub blast_genes{
	my ($hashofseqs,$genemodel)=@_;
		my @hsps=();
		
		
		create_database(\%$hashofseqs);
		
		print "Gene model searching using $genemodel \n";
			
			my @hits=blast_gene_model(\$genemodel);
			push (@hsps,@hits);
		
		
		return @hsps;
	}
	
sub blast_gene_model{
	
	my ($genemodel)=@_;
	print "Using this path for the genemodel sequences ". $$genemodel."\n";
	
	my @params = (program  => 'blastn', database => "./tmp_blast/db" );
	my	$blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
	my $inseq = Bio::SeqIO->new(-file => "<$$genemodel", -format => 'fasta' );
    my @aohits=();
       while (my $seq = $inseq->next_seq) {
       			print $seq->id()."\n";
       	 		my $query_length= $seq->length();
   				#print "Query Length $query_length\n";
       			my  $report_obj = $blast_obj->blastall($seq);
       			my %hits=report_filter($report_obj,$seq);        					
       			#print keys %hits;
       			
       			if (scalar(keys %hits)>0){
       				#print ("RETURNING VALID HITS ", %hits,"\n");
       				push (@aohits,\%hits)
       			}
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
 			 			 while( my $hit = $result->next_hit ) {
   							while( my $hsp = $hit->next_hsp ) {
          								if( $hsp->length('total') >$query_length*0.6 && $hsp->frac_identical( ['query'|'hit'|'total']) >0.99 ) {
     												$hashofhits{$hit->name()}=$seq;
     												#print "Significant hit on scaffold ".;
     												
          								}
    				       			}
			 			 		}
       			 			}

		
		return %hashofhits;
		
}




#for blasting to a reference to annotate
sub blast_orthologs{
	
	my ($inhash)=@_;
	
	my @params = (program  => 'tblastx', database => "home/harrir/projects/gene_identification/arab_blast/db" );
	my $blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
   # my @aohits=();
    my %hash=%$inhash;
    
    my %ho_hits;
      
       foreach (keys %hash){
       			print $_."\n";
       			#exit;
       			my $seq= $hash{$_};
       			#print "Seq ".$seq->length();
       			#exit;
             	my $query_length= $seq->length();
             	my $name=$seq->id();
   				print "Query Length $query_length\n";
       			my  $report_obj = $blast_obj->blastall($seq);
       			my %hits=report_filter_ortho($report_obj,$seq);        					
       			#print keys %hits;
       			
       			if (scalar(keys %hits)>0){
       				#print ("RETURNING VALID HITS ", %hits,"\n");
       				$ho_hits{$_}=\%hits;
       			}
      	}
       	
       	
       	return %ho_hits;
	
}

sub report_filter_ortho{
	my ($report_obj,$seq)=@_;
	my $query_length=$seq->length();
	my $query_name=$seq->id();
	my %hashofhits=();


#	print ("Query Length ", $query_length,"\n");
		while(my $result = $report_obj->next_result ) {
       				my $hitcount= $result->num_hits;
 			 			 while( my $hit = $result->next_hit ) {
 			 			 	my $name=$hit->name();
 			 			 	
   								while( my $hsp = $hit->next_hsp ) {
   									
   		
   										
          								if( $hsp->length('total') >50 && $hsp->frac_identical( ['query'|'hit'|'total']) >0.25 ) {
          										my @arrayofdata=();
          											#print "BLAST REPORT\n";
          											#print "Scaffold $name \t";
          										#	print "Query $query_name \t";
          											push (@arrayofdata,$query_name);
          											push(@arrayofdata,$name);
          											my $hit=$hsp->hit;
          											#print ( "Descriptor ", $hit->seqdesc, "\n");
          											push(@arrayofdata, $hit->seqdesc);
          											#print $hsp->hit_features."\n";
          											#print $hit->annotation()->display_text();
          											#print ("Frac ident ", $hsp->frac_identical( ['query'|'hit'|'total']),"\n" );
          											push(@arrayofdata,$hsp->frac_identical( ['query'|'hit'|'total']));
   												#	print ( "Hit length ", $hsp->length, "\n");
   													push(@arrayofdata,$hsp->length);
   	 											#	print ( "Hit rank ", $hsp->rank,"\n");
   	 												push(@arrayofdata,$hsp->rank);
    											#	print ("Hit start ",$hsp->start('hit'),"\n");
    												push(@arrayofdata,$hsp->start('hit'));
 											     #	print ("Hit end " , $hsp->end('hit'), "\n");
 											     	push(@arrayofdata,$hsp->end('hit'));
 											     #  print ("Strand " , $hsp->strand('hit'), "\n");
 											     	push(@arrayofdata,$hsp->strand('hit'));
     											#	print $hsp->seq( 'hit' )->seq."\n\n";
     												$hashofhits{$name}=\@arrayofdata;
          								}
    				       			}
			 			 		}
       			 			}

		
		return %hashofhits;
}

sub arab_blaster{
	my ($hash,$arab_path)=@_;
	my %inthash=%{$hash};
	
	print "Blasting the hits against arabidopsis\n";
	#print "$$arab_path\n";
	my @argv=("./arab_blast.sh","$$arab_path");
	system(@argv);
	print "BLAST\n";
	#now the blast database is made blast the identified things agaisnt the adab things. 	
	my %ho_hits=blast_orthologs(\%inthash);
	#print keys %ho_hits;
	return %ho_hits;
	
	
}


sub blast_output{
	
	my ($hits,$outputpath,$scaf)=@_;
	print "Printing Summary of Results to $outputpath\n";

open(OUT, ">$outputpath"); #open for write





my %hits=%{$hits};

print OUT ("scaffold\t","fragaria_gene\t","desc_1\t","desc_2\t","frac_ident\t","length\t","rank\t","start\t","end\t","strand\n");

	foreach my $scaffold (keys %hits){
		my %hoh=%{$hits{$scaffold}};
		
		#array of hashes
			foreach  (keys %hoh){
				my @hit=@{$hoh{$_}};
					#print $scaf."\t";	
					print OUT $scaf."\t";
						foreach (@hit){
						#	print $_."\t";
							print OUT $_."\t";
						}	
					#	print "\n";
						print OUT "\n";
						
					}
			}	
	
	close OUT;
	}



sub output_genes{
	
	my ($hits,$outputpath)=@_;
	print "Printing Summary of genes to $outputpath\n";
	open(OUT, ">$outputpath"); #open for write
		
		
		foreach my $hit (@$hits){
			my %hash= %{$hit};
			
			
			foreach (keys %hash){
				print OUT $_."\t";
				print OUT $hash{$_}->id()."\n";
			}
		}
	close OUT;
	}





