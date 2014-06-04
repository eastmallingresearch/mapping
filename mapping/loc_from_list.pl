#!/usr/bin/perl


use warnings;
use strict;

#linkage map (shitty windows control chars)
my $map=shift;
#loc file 
my $loc=shift;
#out file
my $output=shift;
#simple list of extra loci to add to loc file
my $extramap=shift;
#the name of the lg
my $lgname=shift;

#loc_from_list.pl lg3b.txt all.loc lg3b_list.loc extra.txt lg3b



print "TAKING CHROMOSOMAL MAP FILES AND MAKING MASTER LOC FILE \n";
#print "Opening names file $names\n";
my @loc_data= read_table(\$loc); 		
my @parsed_loc=parse_loc(\@loc_data);
my %loc_hash=make_hash(\@parsed_loc);
my @map=read_map(\$map);
my @parsed_map=parse_map(\@map);
my @processed_list=();
my @extramap=read_table(\$extramap);
my @merged = (@parsed_map, @extramap);
my $count=scalar(@merged);

print "COUNT $count \n";

foreach (@merged){
	chomp $_;
	#print $_."\t";
	#looks like there are case issues between the files- needs attention, should match 
	#case independent then use the matched one (i.e. the one from the loc file??
#	print $loc_hash{$_}."\n";
	if (exists $loc_hash{$_}){
		push (@processed_list,$loc_hash{$_});
	}
	else{
	print "Incongruence between master loc and map loc $_ - trying lower case\t";
		$_=~tr/[A-Z]/[a-z]/;
		print $_."\n";
		#exit;
		push (@processed_list,$loc_hash{$_});
		#print $loc_hash{$_};
		
		}
	
}

#exit;

#print @processed_list;
output_marker(\@processed_list, \$output,\$count,\$lgname);

exit;


sub make_hash{
    my ($file)=@_;
    
    my %hashofalleles=();
    
    foreach (@$file){
      	my @split=split("\t",$_);
    	$hashofalleles{$split[0]}=$_;
    }
    
    return %hashofalleles;
    
}


sub read_table{
    my ($file)=@_;
    
    print "Reading in $$file \n ";
    open(IN,$$file);
    my @array= <IN>;
    close IN;
    return @array;
}

sub read_map{
    my ($file)=@_;
    
    print "Reading in $$file \n ";
    open(IN,$$file);
    my $scalar = <IN>;
    close IN;
    $scalar=~tr/\cM/\n/;
    my @array=split("\n",$scalar);
    return @array;
}

sub parse_loc{
my($multiplex_data)=@_;
my $start=7;
my @loc=();
for (my $i=$start;$i<scalar(@$multiplex_data);$i++){
	
	
	if (@$multiplex_data[$i] =~/(\d ; \d)/){
		print @$multiplex_data[$i];
		next;
		}
		elsif (@$multiplex_data[$i] =~/individual names/){
			next;
		}
		else{
	push (@loc, @$multiplex_data[$i]);
		}

}
return @loc;

}

sub parse_map{
my($map_data)=@_;
my $start=1;
my @loc=();
for (my $i=$start;$i<scalar(@$map_data);$i++){
	
	#print $$map_data[$i];
	my @split=split("\t",$$map_data[$i]);
	push (@loc, $split[2]);
}
return @loc;

}

sub output_marker{

my ($combined,$name,$num,$lg)=@_;
my @combined_marker=@{$combined};

print ("NUMBER OF LOCI FOR OUTPUT ", scalar(@combined_marker),"\n");

  open(OUT,">$$name");
  print OUT "name = $$lg\r\n";
  print OUT "popt = CP\r\n";
  print OUT "nind = 188\r\n";
  print OUT "nloc = ".$$num."\r\n";
  print OUT "\r\n";
  
        for (my $i=0; $i<$$num; $i++){
        	#print "OUTPUTING $i \n";
			print OUT $combined_marker[$i]."\n";
		}

print OUT "\r\n";



}














