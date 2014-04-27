#!/usr/bin/perl


use warnings;
use strict;

my $indir=shift;
my $primerdir=shift;
my $seqdir=shift;
my $sumdir=shift;

#./bulk_primer_v2.pl ~/documents/emr/bbsrc_verticillium/vesca/results/candidates/ ~/documents/emr/bbsrc_verticillium/vesca/results/primers/ ~/documents/emr/bbsrc_verticillium/vesca/results/seq/ ~/documents/emr/bbsrc_verticillium/vesca/results/summary/


my @hitlist=<$indir*>;

foreach(@hitlist){
	
	my $hitscaffold=$_;
	$hitscaffold=~s/$indir//;
	$hitscaffold=~s/.txt//;
	print "scaffold being examined ", $hitscaffold,"\n";
	my $newprimerdir=$primerdir.$hitscaffold."/";
	my $newseqdir=$seqdir.$hitscaffold."/";
	my $newsumdir=$sumdir.$hitscaffold."/";
	system("mkdir", $newprimerdir);
	system("mkdir", $newseqdir);
	system("mkdir", $newsumdir);
	#exit;
	
	my @aog=parse_in($_);

	foreach (@aog){
		$_=~s/gene/mrna/;
		my $input=$_.".1";
		my $primer=$newprimerdir.$input."/";
		my $seq=$newseqdir.$input."/";
		my $sum=$newsumdir.$input."/";
		
		
		system("mkdir", $primer);
		system("mkdir", $seq);
		system("mkdir", $sum);
		print "executing ",$input,"\t",$newprimerdir,"\t",$newseqdir,"\n";
		system("./primer_designer_v2.pl", $input ,$primer, $seq, $sum);
		
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



#./primer_designer.pl mrna01837.1 ~/documents/emr/bbsrc_verticillium/vesca/results/primers/  ~/documents/emr/bbsrc_verticillium/vesca/results/seq/


