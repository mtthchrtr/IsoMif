#!/usr/bin/perl

#by Matthieu Chartier 10 juin 2013
#Description
#This program runs a batch of mifView system calls

my $mifViewPath="/Users/matthieuchartier/Development/Isomif/mifView.pl";

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  #print $ARGV[$i]."\n";
  if($ARGV[$i] eq "-e"){ $ex=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-f"){ $dirIn=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-o"){ $dirOut=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeMifView\n##################\n";
    print "-e         <path to mifView.pl>\n";
    print "-f         <dir of input files>\n";
    print "-o         <dir of output files>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

my @mifs=glob $dirIn."/*";
foreach my $file (@mifs){
  next unless($file=~/\.mif$/);
  my $cmd="perl $ex -m $file -o $dirOut";
  print $cmd."\n";
  system($cmd);
}

