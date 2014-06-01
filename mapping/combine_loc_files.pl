#!/usr/bin/perl


use warnings;
use strict;

my $pathtofiles=shift;
my $pathtooutput=shift;
my $names=shift;

print "TAKING CHROMOSOMAL LOC FILES AND MAKING MASTER LOC FILE \n";

print "Opening names file $names\n";

#my $names_file=$pathtofiles.$name1.$i."/";
	#print $d_file."\n";
	opendir DIR, $pathtofiles;
	my @namesfiles = readdir(DIR);
#print @namesfiles;


#./scripts/combine_loc_files.pl ./chromosomal_loc_files/ ./master_loc_files/ ./chromosomal_list/names.txt 

my @total_loc=();
my $count=0;

for (my $i=0;$i<scalar(@namesfiles);$i++){
	
	my $d_file=$pathtofiles.$namesfiles[$i];
	print $d_file."\n";
#	opendir DIR, $d_file;
#	my @files = readdir(DIR);
	
    #	for (my $j=0; $j<scalar(@files); $j++){
  #  		my $p1= $d_file.$files[$j]."\n";
    		my @loc_data= read_table(\$d_file); 		
    		my @parsed_loc=parse_loc(\@loc_data);	    		
    	 	my $add= scalar@parsed_loc;
    		$count=$count + $add;
    		push(@total_loc,\@parsed_loc);
    		
    #	}    	
    }




	print "OUTPUTTING MARKER for $count markers\n";
	my $outpath=$pathtooutput;
	output_marker(\@total_loc, $outpath, 'output', \$count);
	
	

exit;



sub output_marker{

my ($combined, $pathtooutput, $name,$num)=@_;
my @combined_marker=@{$combined};
my $loc_name=$pathtooutput.$name.'.loc'; 


  open(OUT,">$loc_name");
  print OUT "name = EMxFE\r\n";
  print OUT "popt = CP\r\n";
  print OUT "nind = 188\r\n";
  print OUT "nloc = ".$$num."\r\n";
  print OUT "\r\n";
  
        for (my $i=0; $i<scalar(@combined_marker); $i++){
	my @array =$combined_marker[$i];
	
		for(my $j=0;$j<scalar(@array); $j++){
	 	#print $array[$j]."\t";
	    	my @deref = @{$array[$j]};
		#print @deref;
		
			foreach (@deref){
			#print $_;
			print OUT $_."\t";
			}
		
		}
	
	#print OUT $array[scalar(@array)-1];
	
	#print "\n";
	}

print OUT "\r\n";



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



sub read_table{
    my ($file)=@_;
    
    print "Reading in $$file \n ";
    open(IN,$$file);
    my @array= <IN>;
    close IN;
    return @array;
}




