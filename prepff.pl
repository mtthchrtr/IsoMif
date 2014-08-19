#!/usr/bin/perl
use strict;
use warnings;

# my $hivePath="/home/matthieu/hive/";
# my $rootPath="/home/matthieu/Isomif/";
my $hivePath="/Users/matthieuchartier/hive/";
my $rootPath="/Users/matthieuchartier/Development/Isomif/";

#Get atSingle info
my %pseudoToID=();
my %IDtoPseudo=();
my @at=();

#Get atomSingle info
my %aaATID=();
my %aaATOMID=();
my @atom=();
open IN, "<".$rootPath."forcefield_files/ref/atomSingle.isomif" or die "cant open file ./atSingleInfo";
open OUT, ">".$rootPath."forcefield_files/atomSingle.isomif" or die "cant open file ./atSingleInfo";
while(my $line=<IN>){
	print OUT $line;
	unless($line=~/^#/){
		my @info=split(/\s+/,$line);
		$atom[$info[0]][0]=$info[1]; #res
		$atom[$info[0]][1]=$info[2]; #atom
		$atom[$info[0]][2]=$info[3]; #atomtype
		$aaATID{$info[1]}{$info[2]}=$info[3];
		$aaATOMID{$info[1]}{$info[2]}=$info[0];
	}
}
close OUT;
close IN;

#Print getatomid functions header
my @used=();
open OUT, ">".$rootPath."forcefield_files/getAtomId_FAINT.h";
print OUT "int getAtomId_FAINT(string, string);\n\n";
print OUT "int getAtomId_FAINT(string rnam, string anam){\n";
while((my $key)= each %aaATID){
	print OUT "if(rnam.compare(\"$key\")==0){\n";
	while((my $a)= each %{$aaATID{$key}}){
		print OUT "\t if(anam.compare(\"$a\")==0) return $aaATOMID{$key}{$a};\n";
		push @used, $aaATID{$key}{$a} unless($aaATID{$key}{$a}~~@used);
	}
	print OUT "}\n\n";
}
print OUT "\nreturn(-1);}\n";
close OUT;

#Read at single
open IN, "<".$rootPath."forcefield_files/ref/atSingle.isomif";
open OUT, ">".$rootPath."forcefield_files/atSingle.isomif";
printf OUT "%2s %50s %8s %5s %5s %5s %5s %5s", "#Id","Name","Pseudo","Hbdon","Hbacc","Aroma","Chrgd","Hydro";
while(my $line=<IN>){
	next if($line=~/^#/);
	next if($line=~/^$/);
	$line=~s/^\s+//g;
	$line=~s/\s+$//g;
	my @info=split(/\s+/,$line);
	$at[$info[0]][0]=$info[1]; #name
	$at[$info[0]][1]=$info[2]; #pseudo
	$at[$info[0]][2]=$info[3]; #Hb Donnor
	$at[$info[0]][3]=$info[4]; #Hb Acceptor
	$at[$info[0]][4]=$info[5]; #Aromatic
	$at[$info[0]][5]=$info[6]; #Charged amino-acid
	$at[$info[0]][6]=$info[7]; #hydrophobic
	$pseudoToID{$info[2]}=$info[0];
	$IDtoPseudo{$info[0]}=$info[2];
	printf OUT "\n%3d %50s %8s %5d %5d %5d %5d %5d",$info[0],$info[1],$info[2],$info[3],$info[4],$info[5],$info[6],$info[7];
}
close IN;
close OUT;

open IN, "<".$rootPath."forcefield_files/ref/atomSingle.isomif" or die "cant open file ./atSingleInfo";
open OUT, ">".$rootPath."forcefield_files/atoms";
# printf OUT "%10s %10s %10s %10s %10s\n","res","atom","atomType","hBondAngle","armAngle";
while(my $line=<IN>){
	next if($line=~/^$/);
	my @info=split(/\s+/,$line);
	printf OUT "%10s %10s %10s\n",$info[1],$info[2],$at[$info[3]][1];

}
close OUT;
close IN;

#Set atPairs to 0
my @atPairs=();
for(my $i=0; $i<@at; $i++){
	for(my $j=0; $j<=$i; $j++){
		$atPairs[$i][$j]=0;
		$atPairs[$j][$i]=0;
	}	
}

#Read epsilon matrix
open IN, "<".$rootPath."forcefield_files/ref/epsilonMatrix" or die "cant open";
my $count=0;
my @x=();
while(my $line=<IN>){
	chomp($line);
	$line=~s/^\s+//g;
	$line=~s/\s+$//g;
	if($count==0){
		my @info=split(/\s+/,$line);
		for(my $i=0; $i<@info; $i++){
			$x[$i]=$info[$i];
		}
	}else{
		my @info=split(/\s+/,$line);
		my $y=$info[0];
		for(my $i=1; $i<@info; $i++){
			$atPairs[$pseudoToID{$x[$i-1]}][$pseudoToID{$y}]=$info[$i];
			$atPairs[$pseudoToID{$y}][$pseudoToID{$x[$i-1]}]=$info[$i];
		}
	}
	$count++;
}
close IN;

open OUT, ">".$rootPath."forcefield_files/atPair.isomif" or die "cant open atpair";
for(my $i=0; $i<@at; $i++){
	for(my $j=0; $j<=$i; $j++){
		print OUT "$i $j $atPairs[$i][$j]\n";
	}	
}
close OUT;







