#!/usr/bin/perl


use warnings;
use strict;

my $indir=shift;
my $out=shift;

open OUT, ">$out";

#./candidate_blast.pl ~/documents/emr/bbsrc_verticillium/vesca/results/candidates/ ~/documents/emr/bbsrc_verticillium/vesca/results/gs_map/map_candidates.txt


my @hitlist=<$indir*>;

foreach(@hitlist){
	
	my $hitscaffold=$_;
	$hitscaffold=~s/$indir//;
	$hitscaffold=~s/.txt//;
	
	my @aog=parse_in($_);

	foreach (@aog){
		print $hitscaffold,"\t";
		print $_."\n";
		print OUT $hitscaffold,"\t";
		print OUT $_."\n";
		
		}

	}




###subs

sub parse_in {
	my ($in)=@_;

	open IN ,"<$_";
	my @aof = <IN>;
	my @aog=();
	foreach (@aof){
	chomp;
	push (@aog,$_);
	}
	return (@aog);

}


