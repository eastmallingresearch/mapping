#!/usr/bin/perl


use warnings;
use strict;

my $pathtodata=shift;
my $pathtolist=shift;
my $pathtooutput=shift;
my $pathtoloc=shift;
my $name=shift;
my @list=($pathtolist);

#./scripts/format_col.pl ./edited_tables/checked_csv/ ./multiplex_names/processed/ ./transposed_tables/ ./single_loc_files/ multiplex_1

print "Analysing $name \n";
my $fn=$pathtolist.$name.".txt";
my @multiplex_table= read_table($fn);

foreach(@multiplex_table){
    chomp $_;
    print "Retrieving data for ". $_."\n";
    my $fn1=$pathtodata.$_.".csv";
    print "Filename is $fn1\n";
    my @data_table=read_table($fn1);
    #copy last 4 to front and reduce to two

#print @data_table;
#exit;
    my @split_data=split_data(\@data_table);
    my @parents_data=find_parents(\@split_data);
    #re name the first column
    #insert a segregation column
    my @transformed_data=transform_table(\@parents_data);
#print @transformed_data;    
 #foreach (@transformed_data){
	#my @array =@{$_};
	#foreach (@array){
	 #print $_."\t";
   #  print OUT $_."\t";
	#}
 	#print OUT "\n";
	#print "\n";
	#}

#exit;
    my @named_markers=name_markers(\@transformed_data,$_);
   # exit;
   my $primer_name=$pathtooutput.$name."/".$_.'.txt';   
   print "primer_name $primer_name \n";
 # exit;
     open(OUT,">$primer_name");
        foreach (@named_markers){
	my @array =@{$_};
	foreach (@array){
#	 print $_."\t";
     print OUT $_."\t";
	}
 	print OUT "\n";
	#print "\n";
	}
#exit;

my $loc_name=$pathtoloc.$name."/".$_.'.loc';   
my $locs=scalar(@named_markers)-1;
  open(OUT,">$loc_name");
  print OUT "name = EMxFE\r\n";
  print OUT "popt = CP\r\n";
  print OUT "nind = 188\r\n";
  print OUT "nloc = ".$locs."\r\n";
  print OUT "\r\n";
  
        for (my $i=1; $i<scalar(@named_markers); $i++){
	my @array =@{$named_markers[$i]};
	
		for(my $j=2;$j<scalar(@array)-1; $j++){
	 	#print $array[$j]."\t";
	    	print OUT $array[$j]."\t";
		}
	
	print OUT $array[scalar(@array)-1];
	
	print OUT "\r\n";
#	print "\n";
	}



}















sub split_data{
    my ($table)=@_;
    my @AoA=();

    foreach (@$table){
    	#print $_;
    	my $end=substr $_,-1;
    	my $end2=substr $_,-3;
		if($end =~m/,/){
		 substr($_, -1) = ",-"; 
		}
		if($end2 =~m/,\r/){
		 substr($_, -3) = ",-"; 
		}
		
	$_=~s/ //g;
	$_=~s/\n//g;
	$_=~s/\r//g;
	$_=~s/(^|,)(?=,|$)/${1}-/xg;
	#print $_."\n";
	my @split=split(',',$_);
	push(@AoA,\@split);
	}
	#exit;
    return @AoA;
}



sub transform_table{
    my ($AoA)=@_;
    my $i= 0;
    my $j= 0;
    my @transposed=();


for my $row (@$AoA) {
  for my $column (0 .. $#{$row}) {
  #	print $row->[$column]."\n";
    push(@{$transposed[$column]}, $row->[$column]);
   }
 }
 return @transposed;
}


sub read_table{
    my ($file)=@_;
    
    print "Reading in $fn \n ";
    open(IN,$file);
    my @array= <IN>;
    close IN;
    return @array;
}






#this subroutine moves the parents into the correct position
sub find_parents{

my ($data)=@_;

my @deref=@$data;
my @arranged = ();
my $start= scalar(@$data)-4;
my $end= scalar(@$data);
my $length=();
for (my $i=0; $i<$start; $i=$i+1){
	
	if ($i>0){
	my @tmp=@{$deref[$i]};
	$tmp[0]=$i;
	push(@arranged,\@tmp);
	}
	else{
	$length=scalar(@{$deref[$i]});
	push(@arranged,@$data[$i]);
	}
}

my @seg=();
push (@seg,'Seg');

for (my $i=0;$i<$length-1;$i++){
push (@seg,'-');
}

unshift(@arranged,\@seg); 

for (my $i=$start; $i<$end; $i=$i+2){
	unshift(@arranged,@$data[$i]);
}

return (@arranged);

}




#this subroutine names the cluster of markers
sub name_markers{
	my ($trans_data,$pn)=@_;
	my @named_markers=();
	my @deref=@$trans_data;

	#num of markers
	my $j= scalar @{$trans_data};
	push (@named_markers,$deref[0]);
	
	for (my $i=1;$i<$j;$i++){
	
		
		my @marker=@{$deref[$i]};
		my @return_marker=name_marker_type(\@marker,$pn);
		push (@named_markers,\@return_marker);
		

	}

return @named_markers;

}


#this subroutines names coordinates the naming of as single  marker type
sub name_marker_type{
	my ($markers,$pn)=@_;
	#data start
	my $k=4;
	#data end
	my $l=scalar(@{$markers});
	
	my $p1=@$markers[0];
	my $p2=@$markers[1];
	
	my $id=marker_id($p1,$p2);
	@$markers[3]=$id;
	my @renamed=rename_marker(\@$markers,$k,$l,$pn);
	
	return @renamed;
	
}

#this subroutine id's the segregation type
sub marker_id{
my($a,$b)=@_;

	if ($a eq $b){
	return ('<hkxhk>');
	}

	if ($a ne $b && $a eq '-'){
	return ('<nnxnp>');
	}

	if ($a ne $b && $b eq '-'){
	return ('<lmxll>');
	}


}

#this subroutine renames the markers and names the allele

sub rename_marker{
my ($array,$start,$end,$pn)=@_;
my @deref=@$array;
#determine seg type and rename marker and coding
print "MARKER TYPE $pn IS $deref[3] $deref[0] $deref[1]\n";
my @setting=();
if ($deref[3] eq '<lmxll>'){
@setting=('lm','ll','a','--');
#$deref[2]=$setting[2].$pn.'-'.$deref[0];
$deref[2]=$pn.'-'.$deref[0];
#print $deref[2];

}
elsif($deref[3] eq '<nnxnp>'){
@setting=('np','nn','b','--');
#$deref[2]=$setting[2].$pn.'-'.$deref[1];
$deref[2]=$pn.'-'.$deref[1];
#print $deref[2];
}
elsif($deref[3] eq '<hkxhk>'){
@setting=('h-','kk','','--');
#$deref[2]=$setting[2].$pn.'-'.$deref[1];
$deref[2]=$pn.'-'.$deref[1];
}

#rename all alleles
for(my $i=$start; $i<$end;$i++){

if ($deref[$i] =~ /^-?\d+\.?\d*$/){
#print "$deref[$i] int\n";
$deref[$i]=$setting[0];
#print $deref[$i];
}
elsif ($deref[$i] =~ /^-?\w+\.?\w*$/){
#print "$deref[$i] int\n";
$deref[$i]=$setting[3];
#print $deref[$i];
}
else{
#print "$deref[$i] notint\n";
$deref[$i]=$setting[1];
#print $deref[$i];
}

}
return @deref;
}






#unused at present



sub check_ranges{

}
sub sort_table{

    my ($table)=@_;
    my @sorted=();
    my %uniques=number_of_uniques(\@$table);
    print "Number of markers is ", scalar keys %uniques, " \n";
    my $num=scalar keys %uniques;
    #my @roundings=();
    
    
    for (my $i=1;$i<scalar(@$table);$i=$i+4){
	#print $val;
	my @tabulated=();
	
	for(my $j=0; $j<$num;$j++){
	    chomp $$table[$i+$j];
	my @line=split('\t',$$table[$i+$j]);
	my $start=$line[0];
	if ($j==0){
	    #print $line[0]."\t";
	    push(@tabulated,$line[0]);
	}
	    if (exists $line[2]){
		#print $line[2]."\t";
		push(@tabulated,$line[2]);
	    }
	    else{
		push(@tabulated,'');
		#print "\t";
	    }
	}
	print "\n";
	push (@sorted,\@tabulated);
    }
    


    return @sorted;
}


sub number_of_uniques {
    my ($num)=@_;
    my @list;
    
    my $count=0;
    foreach (@$num){
	$count=$count+1;
	my @tmp= split('\t',$_);
	if ($count>1){
	push (@list, $tmp[1]);
	}
    }

    my $txt=join('',@list);

    #my %s;
    #print scalar (grep !$s{$_}++, split '', $_[0]);
    my %same = map { $_, 1 } split //, $txt;
    print scalar keys %same, ' unique: [', keys %same, "]\n";
    foreach (keys %same){
	print ($_ ,"\t", $same{$_},"\n");
    }
    exit;

    my @sizes=();

    foreach my $let (keys %same){
	my $count=0;
	my @nums=();

	foreach (@$num){
	    chomp $_;
	    $count=$count+1;
	    my @tmp= split('\t',$_);
	
	    if ($count>1 && $let eq $tmp[1] && exists $tmp[2]){
		my $round= sprintf("%.0f", $tmp[2]);
		push (@nums, $round);
		    }
	   
	}
	my @uniques=unique(\@nums);
	print @uniques;
	if (scalar @uniques >1){
	    print "fook oop\n";
	    exit;
	}
	else{
	    $same{$let}=$uniques[0];
	}
    }
   
    
    return %same;
    
}

sub unique {
    my ($a)=@_;
    my %seen = ();
    my @result = grep { !$seen{$_}++ } @$a;
    return @result;
}
