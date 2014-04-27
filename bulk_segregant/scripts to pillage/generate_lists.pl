#!/usr/bin/perl


use warnings;
use strict;

my $pathtocsvdata=shift;
my $pathtounscorelist=shift;
my $pathtomultiplexlist=shift;
my $pathtoprocessed=shift;
my $start_multiplex=shift;
my $max_no=shift;

#This script goes through the data and removes those that have failed from the multiplexes and those that have not yet been scored from the multiplexes and appends them to the 


#./scripts/generate_lists.pl ./edited_tables/checked_csv/ ./edited_tables/unscorable_list/ ./multiplex_names/raw/ ./multiplex_names/processed/ 27 49

my $unscore=$pathtounscorelist.'multiplexes.tab';
my @unscore_table= read_table($unscore);


for (my $i=$start_multiplex;$i<=$max_no;$i++){
	my $name=$pathtomultiplexlist.'multiplex_'.$i.'.txt';
	
	my $primer_name=$pathtoprocessed."multiplex_".$i.'.txt';   
     	open(MP,">$primer_name");
     
      	my $unproc_name=$pathtoprocessed.'unscored.txt';   
     	open(RS,">>$unproc_name");
     
     	my $fail_name=$pathtoprocessed.'fail.txt';   
     	open(FS,">>$fail_name");
     
	
	my $outname=$pathtomultiplexlist.'multiplex_'.$i.'.txt';
	print "Analysing $name \n";
	my @multiplex_table= read_table($name);
	
		foreach (@multiplex_table){
		chomp $_;
		my $check=check_table(\$_,\@unscore_table);
		print "RETURN CODE $check\n";
		#exit;
		#ret 0 = normal 1=fail  2=unscore 
		
		if($check==0){
		print MP $_."\n";
		print $_."\n";		
}
		elsif($check==1){
		print FS $_."\n";
		}
		elsif($check==2){
		print RS $_."\n";
		}
		
		#exit;
		
		}

}




sub check_table{
my ($val,$table)=@_;
my $return_code=0;	
	foreach (@$table){
	my @split=split('\t',$_);
	#print "SPLITS $split[0] $split[1] $split[2]\n";	
	chomp $split[0];
	chomp $split[1];
	chomp $split[2];
	$split[1]=~s/"//g; 
	$split[2]=~s/"//g; 
	
		if ($split[1]=~m/$$val/){
			print "MATCHED Checking $$val against $split[1] code $split[2]\n" ;
			chomp $split[2];
			$split[2]=~s/\r//;
			
			if ($split[2] eq 'f'){
			print "MATCHED DUE TO FAIL\n";
			$return_code=1;
			
			}
			elsif($split[2] eq 't'){
			print "MATCHED DUE TO UNSCORE \n";
			$return_code=2;
			}
		
		
		}
		#exit;
	}	
return $return_code;
}


sub read_table{
    my ($file)=@_;
    
    #print "Reading in $fn \n ";
    open(IN,$file);
    my @array= <IN>;
    close IN;
    return @array;
}



