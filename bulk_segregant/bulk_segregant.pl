#!/usr/bin/perl


use warnings;
use strict;
use Statistics::ChiSquare;

my $pathtoloc=shift;
my $pathtomap=shift;
my $pathtopheno=shift;
my $parent=shift;

#./bulk_segregant.pl RGxH_all.loc RGxH_all_FINAL.map m7_control.txt lm



my @loc_data= read_table(\$pathtoloc);
my @parsed_loc=parse_loc(\@loc_data);	 
my @pheno= read_table(\$pathtopheno);
my %pheno_hash=make_hash(\@pheno);
my @transformed_pheno=transform_table(\@pheno);
my @map_data= read_table(\$pathtomap);
#exit;


#print @transformed_pheno;
#exit;
my @AoH=();

for (my $i=0; $i<@parsed_loc; $i++ ){
# for (my $i=186; $i<187; $i++ ){

    my @data=split('\t',$parsed_loc[$i]);
    
    #print "Marker type ".$data[1]."  $i \n";
    
    if ($data[1] eq '<lmxll>'){
        # print "Retrieving data for \t". $data[0]."\t$i\t";
    my @marker_number=number_markers(\$parsed_loc[$i],\$data[1]);
    # print @marker_number;
	my %subset=subset_markers(\@marker_number, $transformed_pheno[0]);
    
    # foreach (keys %subset){
         #   print $_."\n"
    #     }
    # print "\n";
    # foreach (keys %pheno_hash){
    # print $_."\n"
    #}
    my @return=check_locus(\%subset,\%pheno_hash);
        #   print @return;
    my %chi_sq=chi_sq(\@return,\$data[0]);
    
        push(@AoH, \%chi_sq);
    }
}

#print @AoH;

exit;

sub chi_sq{
   my ($array,$locus)=@_;
    my %out=();
    #print $$array[0];
    print $$locus."\t";
    print chisquare($$array[0],$$array[1],$$array[2],$$array[3]);
    $out{$$locus}=chisquare($$array[0],$$array[1],$$array[2],$$array[3]);
    print "\n";
    
    return (%out);
    
}



sub make_hash{
    my ($pheno)=@_;
    my %hash=();
    #ß  print "DATA\n";
    foreach (@$pheno){
         my @tmp=split('\t',$_);
        chomp $tmp[1];
        #print $tmp[1];
        #  print $tmp[0]."\t".$tmp[1]."\n";
         $hash{$tmp[0]}=$tmp[1];
      }
           return %hash;
    
}


sub check_locus{
    my ($subset_hash,$pheno)=@_;
    #print "Checking loci segregation\n";
   
    my %ph=%{$pheno};
    my %sh=%{$subset_hash};
    
    my $total_mp=scalar(keys %ph);
    my $total_ma=0;
    my $positive_1=0;
    my $positive_0=0;
    my $negative_0=0;
    my $negative_1=0;
    my @return=();
    
    foreach (keys %ph ){
        #        print "ACCESSION \t $_ \t marker \t ";
        #     print $sh{$_}."\t pheno \t ";
        #     print $ph{$_}."\n";
        
        if ($sh{$_} eq '-'){
            #  print "MISSING MARKER\n";
            next;
        }
        elsif ($ph{$_} ==0 && $sh{$_}==0){
            $negative_0=$negative_0 +1;
            $total_ma=$total_ma+1;
        }
        elsif ($ph{$_} ==0 && $sh{$_}==1){
            $negative_1=$negative_1 +1;
            $total_ma=$total_ma+1;
        }
        elsif ($ph{$_} ==1 && $sh{$_}==0){
            $positive_0=$positive_0 +1;
            $total_ma=$total_ma+1;
        }
        elsif ($ph{$_} ==1 && $sh{$_}==1){
            $positive_1=$positive_1 +1;
            $total_ma=$total_ma+1;
        }
        elsif ($ph{$_} eq '-'){
            $negative_0=$negative_0 +1;
            $total_ma=$total_ma+1;
        }

        
    }
    push (@return, ($positive_0,$positive_1,$negative_0,$negative_1,$total_mp,$total_ma));
    # print ($positive_0,"\t", $positive_1,"\t",$negative_0,"\t",$negative_1,"\t",$total_mp,"\t", $total_ma, "\n");
    return @return;
}



sub transform_table{
    my ($AoA)=@_;
    my $i= 0;
    my $j= 0;
    my @transposed=();
    
    # print "HERE";
    
    for my $val (@$AoA) {
        my @rows=[split('\t',$val)];
        
        for my $row (@rows) {
            for my $column (0 .. $#{$row}) {
                push(@{$transposed[$column]}, $row->[$column]);
                
            }
            
        }

    }
    # print @transposed;
    return @transposed;
}

sub subset_markers{
    
    my($marker_number,$pheno)=@_;
    
    
    my @reduced_pheno=();
    my %hash;
    #print "PHENO @$pheno scalar is scalar(@$pheno)\n";
    
    #print "MARKER NUMBERS @$marker_number \n";
    #print "SCALAR ". scalar(@$marker_number)." \n";
    
    for (my $i=1; $i<=scalar(@$marker_number);$i++){
        #  print $i."\t".@$marker_number[$i-1]."\n";
        
        for (my $j=0; $j<scalar(@$pheno); $j++){
            if (@$pheno[$j]==$i){
                
                #   print "J is $j I is $i marker is @$marker_number[$i]  pheno is @$pheno[$j]\n";
                
                $hash{$i}=@$marker_number[$i-1];
            }
        }
        
    }
    return %hash;

}



sub number_markers{
    
    my($loc_data,$seg)=@_;
    

    #  print "Segregation type $$seg \n";
   
    my @results=();
    if($$seg=~/<lmxll>/){
        #    print "sending to lmxll\n";
        @results=parse_marker(\$$loc_data,'lm');
    }
    elsif($$seg=~/<nnxnp>/){
        #     print "sending to nnxnp\n";
         @results=parse_marker(\$$loc_data,'np');
	}
    elsif($$seg=~/<hkxhk>/){
        #     print "sending to hkxhk\n";
        # @results=efeg(\@sorted_alleles);
	}
    else{
        #      print "ERROR\n";
    }
    return @results;

    
}

sub parse_marker{
    my($marker,$marker_type)=@_;
    
    #print $$marker;
    my @marker_data=split('\t',$$marker);
    my @numbers=();
    #my %hash=();
    
    for (my $i=3;$i<scalar(@marker_data);$i++){
        # print $marker_data[$i];
        
        if ($marker_data[$i] eq $marker_type){
            # print "HERE";
            push(@numbers,'1');
        
        }
            elsif($marker_data[$i] eq '--'){
           push(@numbers,'-');
        }
        else{
           push(@numbers,'0');
        }
       
    }
        
    return @numbers;
    
    #return %hash;
    # exit;
    
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
    
    #    print "Reading in $$file \n ";
    open(IN,$$file);
    my @array= <IN>;
    close IN;
    return @array;
}



