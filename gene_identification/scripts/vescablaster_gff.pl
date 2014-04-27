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
  my $db     = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                   				 -dsn     => 'dbi:mysql:vescagff', -user => 'harrir');

my $name=shift;
my $out_path=shift;
print $out_path;

open OUT, ">$out_path";

#./vescablaster_gff.pl ~/documents/emr/bbsrc_verticillium/vesca/data/scaffolds.txt ~/documents/emr/bbsrc_verticillium/vesca/results/gs_map/map_gff.txt
#./vescablaster_gff.pl /home/harrir/projects/gene_identification/data/scaffolds.txt /home/harrir/projects/gene_identification/results/map_gff.txt

open IN, "<$name";
my @name=<IN>;
close IN;
foreach (@name){

print "SEARCHING $_ \n";
chomp $_;
 
my @features = $db->get_features_by_location($_);

	foreach  (@features){
		if  ($_->type =~/gene/){
		print $_->seq_id,"\t",$_->name."\n";
		print OUT $_->seq_id,"\t",$_->name."\n";
	
		}
	}
	#my @features = $db->get_features_by_location('scf0513192');

}
close OUT;



				
sub parse_gff{

my ($features)=@_;

my $seqobj=();

 	foreach (@$features){
 					
 					my @feat=();
 					$seqobj = Bio::Seq->new( -display_id => $name, -seq=> $_->seq->seq); 
 					
 					my  $new=Bio::SeqFeature::Generic->new(
 						-primary_tag=>'gene',
                            			-seq_id     => $name,
                            			-start      => $_->start,
                            			-end        => $_->end,
                            			-strand     => $_->strand,
                            			-phase     => 'undef'
                            			);
                            		
                            		push (@feat,$new);      
                            		                         		
                            		my @CDS_only  = $_->get_SeqFeatures('CDS');
					
					foreach (@CDS_only){
											
						my  $new=Bio::SeqFeature::Generic->new(
						-primary_tag=>'CDS',
                            			-seq_id     => $name,
                            			-start      => $_->start,
                            			-end        => $_->end,
                            			-strand     => $_->strand,
                            			-phase     => $_->phase);
                            			
						push (@feat,$new);
					}
					
					 $seqobj->add_SeqFeature(@feat);
				
					
					
					
					
 	};
 
 return $seqobj;
 }
 	
