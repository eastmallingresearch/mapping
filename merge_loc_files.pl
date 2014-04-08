#!/usr/bin/perl


use warnings;
use strict;

my $pathtosingle=shift;
my $pathtorecode=shift;
my $pathtooutput=shift;
my $name=shift;
my $start=shift;
my $end=shift;

#./scripts/merge_loc_files.pl ./single_loc_files/ ./recoded_loc_files/ ./mixed_loc_files/ multiplex_ 1 26
my @total_loc=();
my $count=0;

for (my $i=1;$i<=$end;$i++){
	
	my $d_file=$pathtosingle.$name.$i."/";
	#print $d_file."\n";
	opendir DIR, $d_file;
	my @files = readdir(DIR);
	
    	for (my $j=0; $j<scalar(@files); $j++){
    		my $p1= $d_file.$files[$j]."\n";
    		my @loc_data= read_table(\$p1);
    		my @parsed_loc=parse_loc(\@loc_data);
    		my $add= scalar@parsed_loc;
    		$count=$count + $add;
    		push(@total_loc,\@parsed_loc);
    		
    	}    	
    }


for (my $i=1;$i<=$end;$i++){
	
	my $d_file=$pathtorecode.$name.$i."/";
	#print $d_file."\n";
	opendir DIR, $d_file;
	my @files = readdir(DIR);
	
    	for (my $j=0; $j<scalar(@files); $j++){
    		my $p1= $d_file.$files[$j]."\n";
    		my @loc_data= read_table(\$p1);
    		my @parsed_loc=parse_loc(\@loc_data);
    		my $add= scalar@parsed_loc;
    		$count=$count + $add;
    		push(@total_loc,\@parsed_loc);
    		
    	}    	
    }







	print "OUTPUTTING MARKER for $count markers\n";
	my $outpath=$pathtooutput;
	output_marker(\@total_loc, $outpath, 'output', \$count);
	
	





sub output_marker{

my ($combined, $pathtooutput, $name,$num)=@_;
my @combined_marker=@{$combined};
my $loc_name=$pathtooutput.$name.'.loc'; 
print "inside out\n";

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






#efeg recode

sub efeg{
my ($sorted_alleles)=@_;
my @recode=();

print "Sorting out types\n";
#if it is lmxll-1 or 2 etc

my %code1=("lm"=>"f","h-"=>"","np"=>"g","--"=>'-',"kk"=>"","ll"=>"e","nn"=>"e");
my %code2=("lm"=>"f","ll"=>"e","kk"=>"e","h-"=>"g","--"=>'-');
my %code3=("np"=>"f","nn"=>"e","kk"=>"e","h-"=>"g","--"=>'-');

my $length=();
my $lorn=0;
my @AoA=();
	print "Third time sorted ".scalar(@$sorted_alleles)."\n";
#exit;
	foreach (@$sorted_alleles){
	#print $_;
		my @tmp=split ('\t',$_);
			if($tmp[1] eq '<lmxll>'){
				$lorn=1;
			}
		$length=scalar(@tmp);
		push(@AoA,\@tmp);
	}

#exit;

print "SCALAR AOA ".scalar(@AoA)." $length\n";

if (scalar(@AoA) == 2){
	
	if ($lorn==1){
	print "type lm hk code 2\n";

		@recode=codeefeg(\@AoA,\%code2,\$length);
	}
	else{
	print "type hk np code 3\n";

		@recode=codeefeg(\@AoA,\%code3,\$length);
	}
	unshift(@recode,'<efxeg>');
}

elsif(scalar (@AoA) == 3){
	print "type lm hk np code 1\n";
	@recode=codeefeg(\@AoA,\%code1,\$length);
	unshift(@recode,'<efxeg>');
}
#exit;
return (@recode);
}



#coding subroutine

sub codeefeg{

my ($data, $codeno,$length)=@_;
my %code=%{$codeno};
my @recode=();
print "RE CODING MARKER \n";
#exit;
	for (my $i=2;$i<$$length;$i++){
		my $allelecode=();
		#print "NEXT MARKER \n";
			foreach (@{$data}){
				#print "DATA @{$_}[1] \t";
				
					if(@{$_}[1] eq '<lmxll>'){
						chomp @{$_}[$i];
						@{$_}[$i]=~s/\r//g;
						my $key=@{$_}[$i];
						$allelecode=$allelecode.$code{$key};
				#		print @{$_}[$i]."\t".$code{$key}."\n";
					#exit;
					}
					elsif(@{$_}[1] eq '<nnxnp>'){
					chomp @{$_}[$i];
					@{$_}[$i]=~s/\r//g;
					my $key=@{$_}[$i];
					$allelecode=$allelecode.$code{$key};
				#	print @{$_}[$i]."\t".$code{$key}."\n";
					#exit;
					}
					elsif(@{$_}[1] eq '<hkxhk>'){
					chomp @{$_}[$i];
					@{$_}[$i]=~s/\r//g;
					my $key=@{$_}[$i];
				#	print @{$_}[$i]."\t";
				#	print $code{$key}."\n";
					$allelecode=$allelecode.$code{$key};
				#	print @{$_}[$i]."\t".$code{$key}."\n";
					#exit;
					}
		
			}
	#	print "COMBO ".$allelecode."\n";
		if($allelecode eq '---'){
		$allelecode = '--';
		}
		chomp $allelecode;
	push(@recode,$allelecode);

	
	}

return (@recode);
}






#code up an abcd marker 

sub abcd{
my ($sorted_alleles)=@_;
print "Sorting out types\n";
#if it is lmxll-1 or 2 etc

my %code1=("lm1"=>"a","lm2"=>"b","np1"=>"c","np2"=>"d","--1"=>'--',"--2"=>'--',"ll1"=>"","ll2"=>"","nn1"=>"","nn2"=>"");
my %code2=("lm1"=>"a","ll1"=>"b","np1"=>"c","np2"=>"d","--1"=>'',"--2"=>'-',"ll2"=>"","nn1"=>"","nn2"=>"");
my %code3=("lm1","a"=>"lm2"=>"b","np1"=>"c","nn1"=>"","--1"=>'--',"--2"=>'--',"ll1"=>"","ll2"=>"",,"nn2"=>"");

my @AoA=();
	foreach (@$sorted_alleles){
	my @tmp=split ('\t',$_);
	push(@AoA,\@tmp);
	}
my $lm=0;
my $np=0;
my $length=0;
my @order=();
my @recode=();
#my @order=();

foreach(@AoA){
	if(@{$_}[1] eq '<lmxll>'){
	$lm=$lm+1;
	$length=scalar(@{$_});
	#push(@order,'l');
	}
	elsif(@{$_}[1] eq '<nnxnp>'){
	$np=$np+1;
	#push(@order,'n');
	};

}
print "lm $lm np $np \n";

	if($lm ==2 && $np==2){
	print "ABCD\n";
	@recode=codeabcd(\@AoA,\%code1,\$length);
	unshift(@recode,'<abxcd>')
	#print $length;
	}
	elsif($lm ==2 && $np==1){
	print "ABCN\n";
	@recode=codeabcd(\@AoA,\%code3,\$length);
	unshift(@recode,'<abxcd>')
	#print $length;
	}
	elsif($lm ==1 && $np==2){
	print "ANCD\n";
	@recode=codeabcd(\@AoA,\%code2,\$length);
	unshift(@recode,'<abxcd>')
	#print $length;
	}
return (@recode);
}

#coding subroutine

sub codeabcd{
my ($data, $codeno,$length)=@_;

my %code=%{$codeno};
print "RE CODING MARKER \n";
my @recode=();

	for (my $i=2;$i<$$length;$i++){
	my $lm=0;
	my $np=0;
	my $allelecode=();
	print "NEXT MARKER \n";
		foreach (@{$data}){
		print "DATA @{$_}[1] \t";
			if(@{$_}[1] eq '<lmxll>'){
			$lm=$lm+1;
			chomp @{$_}[$i];
			@{$_}[$i]=~s/\r//g;
			my $key=@{$_}[$i].$lm;
			$allelecode=$allelecode.$code{$key};
			print @{$_}[$i].$lm."\t".$code{$key}."\n";
			}
			elsif(@{$_}[1] eq '<nnxnp>'){
			$np=$np+1;
			print "NP VALUE $np \t";
			chomp @{$_}[$i];
			@{$_}[$i]=~s/\r//g;
			my $key=@{$_}[$i].$np;
			$allelecode=$allelecode.$code{$key};
			print @{$_}[$i].$np."\t".$code{$key}."\n";
			};
		
		}
	print "COMBO ".$allelecode."\n";
	if($allelecode eq '---'){
		$allelecode = '--';
	}
	push(@recode,$allelecode);

	}

return (@recode);
}







#sort abcd from efeg
sub combine_markers{
my($loc_data,$alleles,$seg)=@_;

my @allele_no=split(',',$$alleles);
print ("NO OF ALLELES ", scalar(@allele_no),"\n");
my @sorted_alleles=sort_alleles(\@$loc_data,\@allele_no);
print "Segregation type $$seg \n";
print "Sorted alleles again ".scalar(@sorted_alleles)."\n";
my @results=();
if($$seg=~/<abxcd>/){
print "sending to abcd\n";
@results=abcd(\@sorted_alleles);
}
elsif($$seg=~/<efxeg>/){
	print "sending to efeg\n";
	@results=efeg(\@sorted_alleles);
	}
else{
	print "ERROR";
}
return @results;
}




#pull out relevant alleles from marker
sub sort_alleles{
my ($loc_data,$allele_no)=@_;
my @deref=@{$loc_data};
my @matching_alleles=();

	foreach (@$allele_no){
	print "sort alleles sub searching for $_ \n";	
	my $al_no=$_;
	foreach (@deref){
		my @split=split('\t', $_);
		my $allele= $split[0];
		my @name_split=split('-',$allele);
			if ($name_split[1] == $al_no){
			print "matched\n";
			push (@matching_alleles,$_);
			}
	
		}
	}
	print "Matching alleles ".scalar(@matching_alleles)."\n";
return @matching_alleles;
}





sub parse_loc{
my($multiplex_data)=@_;
my $start=5;
my @loc=();
for (my $i=$start;$i<scalar(@$multiplex_data);$i++){

	push (@loc, @$multiplex_data[$i]);


}
return @loc;

}























#unused





















sub split_data{
    my ($table)=@_;
    my @AoA=();

    foreach (@$table){
    	my $end=substr $_,-3;
		if($end =~m/,/){
		 substr($_, -3) = ",-"; 
		}
	$_=~s/(^|,)(?=,|$)/${1}-/xg;
	$_=~s/\n//g;
	$_=~s/\r//g;
	my @split=split(',',$_);
	push(@AoA,\@split);
	}
    return @AoA;
}





sub transform_table{
    my ($AoA)=@_;
    my $i= 0;
    my $j= 0;
    my @transposed=();


for my $row (@$AoA) {
  for my $column (0 .. $#{$row}) {
  	#print $row->[$column]."\n";
    push(@{$transposed[$column]}, $row->[$column]);
   }
 }
 return @transposed;
}





sub read_table{
    my ($file)=@_;
    
    print "Reading in $$file \n ";
    open(IN,$$file);
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
	print "MARKER RENAMING LENGTH\n";
	print scalar(@renamed);
	
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
$deref[2]=$setting[2].$pn.'-'.$deref[0];
}
elsif($deref[3] eq '<nnxnp>'){
@setting=('nn','np','b','--');
$deref[2]=$setting[2].$pn.'-'.$deref[1];
}
elsif($deref[3] eq '<hkxhk>'){
@setting=('h-','kk','','--');
$deref[2]=$setting[2].$pn.'-'.$deref[1];
}

#rename all alleles
for(my $i=$start; $i<$end;$i++){

if ($deref[$i] =~ /^-?\d+\.?\d*$/){
#print "$deref[$i] int\n";
$deref[$i]=$setting[0];
}
elsif ($deref[$i] =~ /^-?\w+\.?\w*$/){
#print "$deref[$i] int\n";
$deref[$i]=$setting[3];
}
else{
#print "$deref[$i] notint\n";
$deref[$i]=$setting[1];
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
