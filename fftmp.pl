#!/usr/bin/perl
use strict;
use warnings;

my $rootPath="/Users/matthieuchartier/Development/IsoMif/";
my $hivePath="/Users/matthieuchartier/hive/";

my @idtoat=();
my @atsArr=();
open IN, "<".$rootPath."forcefield_files/FlexAID/atomTypes";
open OUT, ">".$rootPath."forcefield_files/FlexAID/probes";
while(my $line=<IN>){
	next if($line=~/^#/);
	$line=~s/^\s+//;
	$line=~s/\s+$//;
	my @s=split(/\s+/,$line);
	$s[0]++;
	$idtoat[$s[0]]=$s[2];
	push @atsArr, $s[2];
	print OUT "$s[2] 0 8 10000";
	print OUT "\n" unless($s[0]==40);
	print "$s[0] $idtoat[$s[0]]\n";
}
close IN;
close OUT;

my %mat=();
open IN, "<".$rootPath."forcefield_files/FlexAID/MC_st0r5.2_6.dat";
while(my $line=<IN>){
	next if($line=~/^$/);
	$line=~s/\s+//g;
	if($line=~/^([0-9]+)[-]([0-9]+)=([-0-9\.]+)$/){
		print "$1 $2 $3\n";
		$mat{$idtoat[$1]}{$idtoat[$2]}=$3;
		$mat{$idtoat[$2]}{$idtoat[$1]}=$3;
	}
}
close IN;

open OUT, ">".$rootPath."forcefield_files/FlexAID/epsilons";
printf OUT "%8s","";
foreach my $k1 (@atsArr){ printf OUT " %8s",$k1; }
print OUT "\n";
foreach my $k1 (@atsArr){
	printf OUT "%8s",$k1;
	foreach my $k2 (@atsArr){
		if($k1 eq $k2){
			printf OUT " %8d",0;
		}else{
			print $k1." ".$k2."\n";
			if($mat{$k1}{$k2}==0){
				printf OUT " %8d",0;
			}else{
				printf OUT " %8.3f",$mat{$k1}{$k2};				
			}
		}
	}
	print OUT "\n" unless($k1 eq "solvent");
}
close OUT;















