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

#primer3_core wouldn't recognize tab of SEQUENCE but SEQUENCE TEMPLATE. Solution, download an earlier version of primer3_core (1.1.4 version), make it and change the path 

my $db     = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                   				 -dsn     => 'dbi:mysql:vescagff', -user => 'harrir');
my $name=shift;
my $out_path=shift;
my $seq_out=shift;
my $summary_out=shift;


#./primer_designer_v2.pl mrna06522.1 ~/documents/emr/bbsrc_verticillium/vesca/results/primers/  ~/documents/emr/bbsrc_verticillium/vesca/results/seq/ ~/documents/emr/bbsrc_verticillium/vesca/results/summary/

print "Retrieving feature $name \n";
my $summary=$summary_out.$name.".txt";
print "Summary output path $summary \n";
open SUM, (">$summary");


my @features = $db->get_features_by_name($name);
my $seq=parse_gff(\@features,\*SUM);
#scan for ssrs and add them as objects to the sequence
print "SCANNING FOR SSRs\n";
ssr_scan($seq,\*SUM);
#sequence features of interest for primer design
my @locations=region_locations(\$seq,\*SUM);
#get some coordinates and set max size- for introns- primers are in coding seq
my @coordinates=pick_coordinates(\@locations,'400',$seq);
#exit;
# design the primers and output the primer and selected location sequences
my $coordset=0;
	foreach(@coordinates){
		$coordset=$coordset+1;
		my @pair=@{$_};
		print "RAW COORDINATES " ,$pair[0],"\t",$pair[1],"\t",$pair[2],"..",$pair[3],"\n";
			if ($seq->length<$pair[1]){
				print "TRUNCATION DUE TO LENGTH UNDER/OVERSPILL\n";
				$pair[1]= $seq->length;
				}
			if($pair[0]<1){
				$pair[0]=1;
				}
		print "TRUNCATING SEQUENCE FOR PRIMER DESIGN $pair[0] $pair[1] excluding  $pair[2] $pair[3] \n";
		print SUM "\nSEQUENCE FOR PRIMER DESIGN $pair[0] - $pair[1] \t EXCLUDING $pair[2] - $pair[3] \n";
		#exit;
		my $subseq=$seq->trunc($pair[0],$pair[1]);
		my $results=design_primer($subseq,($pair[2]-$pair[0]),($pair[3]-$pair[0]));
		output_results($results,$name,$out_path,$coordset,$seq,\*SUM);
	}
	
# output the master 'annotated' sequence
output_seq($seq,$seq_out);
close (SUM);

#####SUBROUTINES


#pick primer coordinates

sub pick_coordinates {

my ($loc,$max,$seq)=@_;
print "retrieving location information of primer coords\n";
my @transcript= $seq->get_SeqFeatures();
my @introns=$transcript[0]->introns;

my @all_coord=();
	foreach (@$loc){
		my @tmp=picker(\@{$_},$max,\@introns);
		if (@tmp){
		push(@all_coord,\@tmp);		
		}
	}

return (@all_coord);
}

#do the actual picking
sub picker{
use Math::Round;

my ($loc,$max,$introns)=@_;
my @tmp=@{$loc};
my @picked=();

	
	#exclude very large introns--- make sure to check that the start is not too far back into another intron.. 
	
	if ($tmp[2]>$max){
		#print "this is too big\n";
		return();
	}

	else{
	my $midpoint=$tmp[2]/2;
	print "FEATURE POS $tmp[0] $tmp[1] and length $tmp[2]\n";
	print "MIDPOINT of feature ",$midpoint,"\n";
	my $sub=round(($max-$midpoint)/2,0);
	#print $sub."\t";
	#print $tmp[0]-$sub."\t";
	
	
	my $prelim_s=$tmp[0]-$sub;
	my $prelim_e=$tmp[1]+$sub;
	print " PRELIM COORDS= $prelim_s $prelim_e \n";
		foreach(@$introns){
		
		if ($prelim_s >$_->start && $prelim_s <$_->end){ 
			print "OVERLAPS ADJECENT INTRON 5' OR IS STILL WITHIN INTRON ", $_->start,"..",$_->end ,"- resetting to ",$_->end +1,"\n";
						
					$prelim_s=$_->end+1;	
						}
			
		if ($prelim_e >$_->start && $prelim_e <$_->end){ 
			print "OVERLAPS ADJECENT INTRON 3' OR IS STILL WITHIN INTRON ",$_->start,"..", $_->end, "- resetting to ",$_->start -1,"\n";
			$prelim_e=$_->start-1;
			}
			
		}
	
	if ($prelim_e<$prelim_s || $prelim_s>$tmp[0]){
	print "Buggered this one up, resetting to original values $tmp[0]-$sub $tmp[1]+$sub!\n";
	$prelim_s=$tmp[0]-$sub;
	$prelim_e=$tmp[1]+$sub;
	
	}
	
	push(@picked, $prelim_s);
	#print $tmp[1]+$sub."\n";
	
	print "PUSHING $prelim_s and $prelim_e \n";
	push(@picked, $prelim_e);
	push(@picked,$tmp[0]);
	push(@picked,$tmp[1]);
	}
return @picked;

}


#work out decent locations for primers

sub region_locations{
my ($seq,$fh)=@_;

my @transcript= $$seq->get_SeqFeatures();
#my @introns = $transcript[0]->introns();
#my $cds = $transcript[0]->cds();
#map {print $_->start,'..',$_->end,"\t",($_->end)-($_->start),"\n"} $transcript[0]->introns;
#my @size=map {$_->end-$_->start} $transcript[0]->introns;

print "INTRON LOCATIONS \n";
print $fh "\nINTRON LOCATIONS (START\tEND\tLENGTH\tCLASS) \n";
my @aof=();
my $intron_marker=0;

	if (defined($transcript[0]->introns) ){
		$intron_marker=1;
	};



if ($intron_marker>0){
	foreach($transcript[0]->introns){
		my @intron=();
		print $_->start,"\t", $_->end,"\n";
		print $fh $_->start,"\t", $_->end,"\t";
		push (@intron,$_->start);
		push (@intron,$_->end);
		push(@intron,($_->end)-($_->start));
		print $fh $_->end- $_->start,"\t";
		push(@intron,'intron');
		print $fh "intron\n";
		push (@aof,\@intron);
	}

}

print "SSR LOCATIONS \n";
print $fh "SSR LOCATIONS (START\tEND\tLENGTH\tCLASS\tNESTING) \n";
	for (my $i=1;$i<@transcript;$i++){
	my @ssr=();
	print $transcript[$i]->start,"\t",$transcript[$i]->end,"\t",$transcript[$i]->end-$transcript[$i]->start,"\n";
	print $fh $transcript[$i]->start,"\t",$transcript[$i]->end,"\t",$transcript[$i]->end-$transcript[$i]->start,"\t";
	push (@ssr,$transcript[$i]->start);
	push (@ssr,$transcript[$i]->end);
	push (@ssr,$transcript[$i]->end-$transcript[$i]->start);
	print $fh "ssr \t";
	push (@ssr,'ssr');
	
		if ($intron_marker ==1){
			foreach($transcript[0]->introns){
				if ($transcript[$i]->start > $_->start && $transcript[$i]->start < $_->end){
					print $fh "within intron \n";
					push (@ssr,'intron');
				}
		
			}
				if ($ssr[4] ne 'intron'){
					$ssr[4]= 'exon';
					print $fh "within exon \n";
				}
		}

		if ($intron_marker ==0){
			$ssr[4]= 'exon';
			print $fh "exon \n";
		}


	push(@aof,\@ssr);

}	

return (@aof);
}


#scan for ssrs
sub ssr_scan{

my ($seq,$fh)=@_;
my $seqtot=$seq->seq;
print $fh "\nSSR SCAN (SSR\tMOTIF\tSTART\tEND\t)\n";
	 while ( $seqtot =~ /(([ATGC]{2,})\2{3,})/gi ) {
	   
	    print $1."\t MOTIF \t";
	    print $2."\n";
	    print $fh $1."\t";
	    print $fh $2."\t";
	    
	    my $pos2 = index($seqtot,$1);
	    print $pos2."\t";
	     print $fh $pos2."\t";
	    print $pos2+length($1)."\n";
	    print $fh $pos2+length($1)."\n";
	    
	    
	     my $feat = Bio::SeqFeature::Generic->new( 
		    -start        => $pos2, 
		    -end          => $pos2+length($1),
		    -primary      => 'ssr', # -primary_tag is a synonym
		    -source_tag   => 'ssrscan',
		    -display_name => $2
		    );
	    $seq->add_SeqFeature($feat);
	}

}


#parse the database accession
			
sub parse_gff{
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;


my ($features,$fh)=@_;
my $transcript=Bio::SeqFeature::Gene::Transcript->new();
my $seqobj=();

#strand and seqene
my $strand=();
my @CDS_only=();

	foreach (@$features){
		$strand= $_->strand;
		$seqobj = Bio::Seq->new( -display_id => $name, -seq=> $_->seq->seq, -strand=>$strand); 
		@CDS_only = $_->get_SeqFeatures('CDS');
	}
	
 print "THIS GENE IS ON THE $strand STRAND\n";
 print $fh "THIS GENE IS ON THE $strand STRAND\n";

my @ordered=();
my $start=();	
# order the exons
	
	if ($strand ==1){
	@ordered= sort { $a->start cmp $b->start} @CDS_only;
	$start=$ordered[0]->start;
	print $start."\n";
	print $ordered[0]->end."\n";
	}
	else{
	@ordered= sort { $b->end cmp $a->end} @CDS_only;
	$start=$ordered[0]->end+1;
	}
	
	
 	foreach (@ordered){
 		
 		my $ne=$start-$_->start;
 		my $ns=$start-$_->end;
 		
 		if ($strand ==1){
 		$ns=$_->start-$start+1;
 		$ne=$_->end-$start+1;
 		}
 		print "Resetting the coordinates\n";
 		
 		print "EXON_START ", $_->start,"\tEND\t", $_->end,"\tNEW_START\t",$ns,"\tNEW END\t",$ne,"\n";
 		print $fh "EXON_START ", $_->start,"\tEND\t", $_->end,"\tNEW_START\t",$ns,"\tNEW END\t",$ne,"\n";
 		
 		my  $new=Bio::SeqFeature::Gene::Exon->new(
							-primary_tag=>'exon',
				  			-seq_id     => $name,
				   			-start      => $ns,
				    			-end        => $ne,
				    			-strand     => $_->strand,
				    			-phase     => $_->phase
				    			);
						       			
							$transcript->add_exon($new,);
				}
		
		 	
	$seqobj->add_SeqFeature($transcript);
	#exit;
	return $seqobj;
 
 
 
 }
 	





#output a basic sequence- this can be improved
sub output_seq{
my ($seq, $path)=@_;

my @transcript= $seq->get_SeqFeatures();
my $cds = $transcript[0]->cds();

my $strand= $transcript[0]->strand;
#print "STRAND $strand\n";

if ($strand==(-1)){
$cds= $cds->revcom();

}

my $id=$seq->id;
my $id2=$id;
$id2=~s/mrna/gdna/;
my $named=$path.$id;
my $named2=$path.$id2;


my $seqout=Bio::SeqIO->new(-file=>">$named2.gbk",
                                -format=>'genbank');
    		 # print the sequence out
    		 	$seqout->write_seq($seq);
    		 	
$seqout=Bio::SeqIO->new(-file=>">$named2.fa",
                                -format=>'fasta');
    		 # print the sequence out
    		 	$seqout->write_seq($seq);
    		 	
$seqout=Bio::SeqIO->new(-file=>">$named.cds.fa",
                                -format=>'fasta');
    		 # print the sequence out
    		 	$seqout->write_seq($cds);

}




#output primers

sub output_results{
my ($results,$name,$out_path,$coordset,$seq,$fh)=@_;

	my $all_results = $results->all_results;
	 #print "ALL the results\n";
 		#foreach my $key (keys %{$all_results}) {
    		#	print "$key\t${$all_results}{$key}\n";
 		#}
	my $count=0;
 	while (my $primed_seq  = $results->next_primer()) {
 		
 		
 		#$seq->add_SeqFeature($primed_seq);
 		$count=$count+1;
 		#print $primed_seq;
 		my $named=$out_path.$name.'_'.$coordset.'_'.$count.".gbk";	
 		my $named1=$out_path.$name.'_'.$coordset.'_'.$count.".fa";
 		my $named_seq=$out_path.$name.'_'.$coordset.'_'.$count.".prm";
 		 		
 		 my $seqout=Bio::SeqIO->new(-file=>">$named",
                                -format=>'genbank');
    		 # print the sequence out
    		 	my $aso = $primed_seq->annotated_sequence();
     			$seqout->write_seq($aso);
     		
     		$seqout=Bio::SeqIO->new(-file=>">$named1",
                                -format=>'fasta');
    			$seqout->write_seq($aso);		
     		
     		
     			my %hash=%{$primed_seq};
     			#print keys %hash;
     			my $tmp= $hash{'annotated_sequence'};
     			#print $tmp->get_SeqFeatures;
     			my @primer_feat=$tmp->get_SeqFeatures;
     			
     			foreach (@primer_feat){
	     			print $fh "START\t".$_->start."\t";
	     			print $fh "END\t".$_->end."\n";
     			
     				print  "START\t".$_->start."\t";
	     			print  "END\t".$_->end."\n";
     			}
     			#exit;
     			open OUT, ">$named_seq";
     			print OUT "LEFT\t",$hash{'left_primer'}->seq->seq;
     			print OUT "RIGHT\t",$hash{'right_primer'}->seq->seq;
     			
     			print  "LEFT\t",$hash{'left_primer'}->seq->seq."\n";
     			print  "RIGHT\t",$hash{'right_primer'}->seq->seq."\n";
     			
     			print $fh "LEFT\t",$hash{'left_primer'}->seq->seq."\n";
     			print $fh "RIGHT\t",$hash{'right_primer'}->seq->seq."\n";
     			close OUT;
     			
     			    			
 	}

}





sub design_primer{
	my ($seq,$exc1,$exc2)=@_;
	  use Bio::Tools::Run::Primer3;
	  use Bio::SeqIO;
	  my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq,
		                                      -outfile => "temp.out",
		                                      -path => "/home/harrir/prog/primer3-1.1.3/");

	print "SEARCHING FOR PRIMERS\n";
	my $len=$exc2-$exc1;
	#print "EXCLUDED REGION $exc1 .. $exc2 length $len\n";
	my $exc=($exc1.",".$len);
	print "EXCLUDED REGION $exc\n";
	
	  # what are the arguments, and what do they mean?
#	  my $args = $primer3->arguments;

	 # print "ARGUMENT\tMEANING\n";
#	  foreach my $key (keys %{$args}) {print "$key\t", $$args{$key}, "\n"}
	#
	  # set the maximum and minimum Tm of the primer
	 	$primer3->add_targets('PRIMER_MIN_TM'=>57, 'PRIMER_MAX_TM'=>63,'PRIMER_PRODUCT_OPT_SIZE'=>300,'PRIMER_OPT_SIZE'=>22,'PRIMER_GC_CLAMP'=>2,'PRIMER_PRODUCT_SIZE_RANGE'=>'200-490','EXCLUDED_REGION'=>$exc);

	  # design the primers. This runs primer3 and returns a 
	  # Bio::Tools::Run::Primer3 object with the results
	 
	 my $results = $primer3->run;

	  # see the Bio::Tools::Run::Primer3 pod for
	  # things that you can get fro22m this. For example:

	  print "There were ", $results->number_of_results, " primers\n";
#exit;
	return $results;

	
}








    					
    
    
    
    
    
    
    		
 	

