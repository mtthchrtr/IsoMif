#!/usr/bin/perl
use strict;
use warnings;

#by Matthieu Chartier
#Description
#This program generates the mifView pml files

my @probesInt=();
my @ca=();
my @gridInt=();
my @probes=();
my @grid=();
my @probesLab=("HYD","ARM","DON","ACC","POS","NEG");
my @pbColors=("aquamarine","brightorange","blue","red","limegreen","lightmagenta");
my @gridLab=("200","150","100","050");
my @gridColors=("br2","palegreen","hydrogen","cyan");
my $mifViewFolder="./";
my $mifName="";
my $mifFolder="";
my $mifFile="";

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  if($ARGV[$i] eq "-m"){ $mifFile=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-o"){ $mifViewFolder=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeIsoMifView\n##################\n";
    print "-m         <path to mif file>\n";
    print "-o         <mifView output directory>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

if($mifFile eq ""){
  print "Error: Missing mif File\n";
  exit;
}

# $mifViewFolder=&get_dirs("/Users/matthieuchartier/hive/","mifView") if($mifViewFolder eq "");

if($mifFile=~/\/([a-z0-9_-]+)\.mif$/i){
  $mifName=$1;
  $mifFolder=$mifFile;
  $mifFolder=~s/\/$mifName\.mif$//;
}elsif($mifFile=~/^([a-z0-9_-]+)\.mif$/i){
  $mifName=$1;
  $mifFolder="./";
}else{
  print "\nError. Mif filename must contain a-z 0-9 _ or - characters.\n";
  exit;
}

open IN, "<".$mifFile or die "cant open mif file $mifFile";
while(my $line=<IN>){
  if($line=~/^#ATOM/){
    if($line=~/#ATOM\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]{1})\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([0-9]{1})\s+([0-9]{1})$/i){
      push @ca, "$1 $2 $3 $4 $5 $6 $7 $8" if($10==1 && $3 eq "CA");
      # print "$1 $2 $3 $4 $5 $6 $7 $8\n";
    }
    next;
  }elsif($line=~/^#/){
    next;
  }elsif($line=~/^$/){
    next;
  }
  $line=~s/^\s+//g;
  $line=~s/\s+$//g;
  my @info=split(/\s+/,$line);

  #Store vrtx potential interaction
  for(my $i=3; $i<9; $i++){
    push @{$probes[$i-3]}, ($info[0],$info[1],$info[2],$info[9],$info[10],$info[11],$info[12]) if($info[$i]==1);
  }

  #Store vrtx grid presence
  for(my $i=9; $i<13; $i++){
      push @{$grid[$i-9]}, ($info[0],$info[1],$info[2]) if($info[$i]==1);
  }
}
close IN;

open NPML, ">".$mifViewFolder."/".$mifName.".pml";
print NPML "feedback disable,all,output\n";

open IN, "<".$mifFolder."/".$mifName."_cpy.pdb";
print NPML "cmd.read_pdbstr(\"\"\"";
while(my $line=<IN>){
  chomp($line);
  print NPML $line."\\\n";
}
print NPML "TER \\\n\"\"\",\"".$mifName."\")\nshow cartoon\nremove (resn HOH)\nshow sticks, HET & ".$mifName."\nset connect_mode,1\n";
close IN;

# open OUT, ">".$mifViewFolder."/".$mifName.".pdb" or die "Cant open mifView file";
my $it=0;
for(my $i=0; $i<6; $i++){ #Loop each probe
  if($#{$probes[$i]}){
    for(my $g=0; $g<4; $g++){ #Loop each grid resolution
      $probesInt[$i][$g][0]=$it;
      print NPML "cmd.read_pdbstr(\"\"\"";
      for(my $j=0; $j<@{$probes[$i]}; $j+=7){ #For each node
        if($probes[$i][$j+3+$g]==1){ #If its in this grid resolution
          # printf OUT "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\n",$it,$probesLab[$i],$probes[$i][$j],$probes[$i][$j+1],$probes[$i][$j+2];          
          printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,$probesLab[$i],$probes[$i][$j],$probes[$i][$j+1],$probes[$i][$j+2];
          $it++;
        }
      }
      print NPML "TER \\\n\"\"\",\"".$probesLab[$i].$gridLab[$g]."\")\n";
      $probesInt[$i][$g][1]=$it-1;
    }
  }
}
for(my $i=0; $i<3; $i++){
  if(@{$grid[$i]}){
    $gridInt[$i][0]=$it;
    print NPML "cmd.read_pdbstr(\"\"\"";
    for(my $j=0; $j<@{$grid[$i]}; $j+=3){
      # printf OUT "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\n",$it,$gridLab[$i],$grid[$i][$j],$grid[$i][$j+1],$grid[$i][$j+2];          
      printf NPML "HETATM%5d  N   %3s A0000    %8.3f%8.3f%8.3f  0.00 10.00           N\\\n",$it,$gridLab[$i],$grid[$i][$j],$grid[$i][$j+1],$grid[$i][$j+2];          
      $it++;
    }
    print NPML "TER \\\n\"\"\",\"".$gridLab[$i]."\")\n";
    $gridInt[$i][1]=$it-1;
  }
}
# close OUT;

print NPML "feedback enable,all,output\norient\n";

# open PML, ">".$mifViewFolder."/".$mifName.".pml" or die "Can't open file";
# system("cp ".$mifFolder."/".$mifName."_cpy.pdb ".$mifViewFolder."/") unless($mifFolder eq $mifViewFolder);
# print PML "load ".$mifViewFolder."/".$mifName."_cpy.pdb\nshow cartoon\nhide lines\n";
# print PML "remove (resn HOH)\nhide everything, hydrogens\n";
# print PML "load ".$mifViewFolder."/".$mifName.".pdb\n";
for(my $i=0; $i<6; $i++){
  if(@{$probes[$i]}){
    for(my $g=0; $g<4; $g++){
      if($probesInt[$i][$g][0]!=$probesInt[$i][$g][1]){
        # print PML "create ".$probesLab[$i].$gridLab[$g].", id $probesInt[$i][$g][0]-$probesInt[$i][$g][1] & ".$mifName."\n";
        # print PML "show spheres, ".$probesLab[$i].$gridLab[$g]."\nset sphere_scale,0.2,".$probesLab[$i].$gridLab[$g]."\nrebuild\n";
        # print PML "color $pbColors[$i],".$probesLab[$i].$gridLab[$g]."\nhide nonbonded,".$probesLab[$i].$gridLab[$g]."\n";
        print NPML "show spheres, ".$probesLab[$i].$gridLab[$g]."\nset sphere_scale,0.2,".$probesLab[$i].$gridLab[$g]."\nrebuild\n";
        print NPML "color $pbColors[$i],".$probesLab[$i].$gridLab[$g]."\nhide nonbonded,".$probesLab[$i].$gridLab[$g]."\n";
      }
    }
  }
}
for(my $i=0; $i<3; $i++){
  if(@{$grid[$i]}){
    if($grid[$i][0]!=$grid[$i][1]){
      # print PML "create ".$gridLab[$i].", id $gridInt[$i][0]-$gridInt[$i][1] & ".$mifName."\n";
      # print PML "color $gridColors[$i],".$gridLab[$i]."\nshow nonbonded,".$gridLab[$i]."\n";
      print NPML "color $gridColors[$i],".$gridLab[$i]."\nshow nonbonded,".$gridLab[$i]."\n";
    }
  }
}

#Print Ca atoms
# my $nbca=@ca;
# if($nbca>0){
#   print PML "create ca, id ";
#   for(my $i=0; $i<@ca; $i++){
#     my @s=split(/\s+/,$ca[$i]);
#     print PML "$s[3]";
#     print PML "+" unless($i==$#ca);
#   }
#   print PML " & ".$mifName."_cpy\n";
# }
# print PML "color deepteal, ca\nset sphere_scale, 0.3, ca\nshow spheres, ca\n";

# print PML "delete ".$mifName."\nhide lines\nshow sticks, HET & ".$mifName."_cpy\n";
# print PML "save ".$mifViewFolder."/".$mifName.".pse\n";
# close PML;

########################################
#   SUBS
########################################

sub get_dirs{
  my $rootDir=$_[0];
  my $base=$_[1];
  my $listDir=$rootDir.$base."/";
  my $dirString;

  #Get the mifs directory
  opendir my $dh, $listDir or die "$0: opendir: $!";
  my @dirs=();
  my $count=0;
  print "\n\n";
  while (defined(my $name = readdir $dh)) {
    next unless -d $listDir."$name";
    next if($name eq "..");
    print "$count $name\n";
    push @dirs, $name;
    $count++;
  }
  print "\n".$listDir;
  print "\nEnter the number of the desired dir:";
  my $answer=<STDIN>;
  if ($dirs[$answer] eq ".") {
    $dirString=$listDir;
    } else {
      $dirString=$listDir.$dirs[$answer]."/";
    }

    print "Path is: $dirString";
    return($dirString);
  }