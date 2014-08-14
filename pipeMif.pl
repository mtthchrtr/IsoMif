#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager; #Call the package to use for forks.

#by Matthieu Chartier 10 juin 2013
#Description
#This program runs a batch of mif system calls on the jobs list with the mifParameters

#Default values
my $mifPath="";
my $mifParamF="";
my $mifJobsF="";
my $cmdMode=0;
my $tag="pipeMif";
my $jobsDir="";
my $batch=30;
my @mifParam=();
my @cases=();
my @cmds=();

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  #print $ARGV[$i]."\n";
  if($ARGV[$i] eq "-e"){ $mifPath=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-p"){ $mifParamF=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-j"){ $mifJobsF=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-c"){ $cmdMode=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-b"){ $batch=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-t"){ $tag=$tag.$ARGV[$i+1]; }
  if($ARGV[$i] eq "-f"){ $jobsDir=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeMif\n##################\n";
    print "-e         <path to MIF program>\n";
    print "-p         <file with all the isomif parameters>\n";
    print "-j         <file with the list of jobs>\n";
    print "-c         <only print command line [nrg (batch qsub on nrg), 1 (print cmds only), 0 (run locally)]>\n";
    print "-b         <Nb of cmds per job file>\n";
    print "-t         <tag of this job batch>\n";
    print "-f         <dir where to print the job files>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

&storeParams();

&storeCases();

&recur(0,"",0);

&runCmds();

sub runCmds{
  if($cmdMode eq "nrg"){
    my $filenb=0;
    my $count=0;
    open OUT, ">".$jobsDir.$filenb.".pbs" or die "cant open".$jobsDir.$filenb.".pbs";
    print OUT "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -N ".$tag.$filenb."\n";
    for(my $i=0; $i<@cmds; $i++){
      # if($count==ceil(@cmds/$batch)){
      if($count==$batch){
        $count=0;
        close OUT;
        $filenb++;
        open OUT, ">".$jobsDir.$filenb.".pbs" or die "cant open".$jobsDir.$filenb.".pbs";
        print OUT "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -N ".$tag.$filenb."\n";
      }
      print OUT "$cmds[$i]";
      $count++;
      print OUT "\n" unless($count==$batch);
    }
    close OUT;

    my @files=glob($jobsDir."/*");
    foreach my $file (@files){
      print "qsub $file\n";
      system("qsub $file");
      # print "qsub $file\n";
      sleep(0.1);
    }
    &areJobsDone("nrg",$tag);
  }elsif($cmdMode eq "print"){
    for(my $i=0; $i<@cmds; $i++){
      print $cmds[$i]."\n";
    }
  }elsif($cmdMode eq "local"){
    my $nProcess = $batch; # Number of threads allowed simultaneously
    my $fork= new Parallel::ForkManager($nProcess); # The objects that manages the forks.
    for(my $i=0; $i<@cmds; $i++){
      $fork->start and next; # Forking a new child process.
      system($cmds[$i]." > /dev/null 2>&1");
      $fork->finish; # do the exit in the child process.
    }
    $fork->wait_all_children; # Wait for all forks to exit.
  }
}

sub storeCases{
  open IN, "<".$mifJobsF or die "Cant open ".$mifJobsF;
  while(my $line=<IN>){
    if($line!~/^$/){
      chomp($line);
      push @cases, $mifPath." ".$line;
    }
  }
  close IN;
}

sub recur{
  my $level=$_[2];
  $level++;
  for(my $p=$_[0]; $p<@mifParam; $p++){
    for(my $i=0; $i<@{$mifParam[$p][1]}; $i++){
      my $cmd=$_[1]." ".$mifParam[$p][0]." ".$mifParam[$p][1][$i]; 
      if($p==$#mifParam){
        if($level==@mifParam){
          foreach my $c (@cases){
            push @cmds, $c.$cmd;
          }
        }
      }else{
        &recur($p+1,$cmd,$level);
      }
    }
  }
}

sub storeParams{
  my $p=0;
  open IN, "<".$mifParamF or die "Can't open ".$mifParamF;
  while(my $line=<IN>){
    next if($line=~/^$/);
    my @sub=split(/\s+/,$line);
    $mifParam[$p][0]=$sub[0];
    for(my $i=1; $i<@sub; $i++){
      $mifParam[$p][1][$i-1]=$sub[$i];
    }
    $p++;
  }
  close IN;
}

sub areJobsDone{
  my $sys=$_[0];
  my $exitLoopR=0;
  my $exitLoopQ=0;
  my $getout=0;
  my $time=0;
  print "\n\nWaiting for jobs to terminate..";
  while(1){
    sleep 5;
    $time+=5;
    my $string;
    if($sys eq "nrg"){
      $string="qstat | egrep '".$_[1]."'| egrep ' R ' | wc -l |";
      open COM, $string or die "cant open qstat grep check";
      while(my $line=<COM>){
        print "Running: $line\n";
        if($line=~/^0$/){
          $exitLoopR=1;
          last;
        }
      }
      close COM;
      $string="qstat | egrep '".$_[1]."'| egrep ' Q ' | wc -l |";
      open COM, $string or die "cant open qstat grep check";
      while(my $line=<COM>){
        print "Queud: $line\n";
        if($line=~/^0$/){
          $exitLoopQ=1;
          last;
        }
      }
      close COM;
      $getout=1 if($exitLoopR==1 && $exitLoopQ==1);
    }

    last if($getout==1);
  }

  print "Job ".$_[0]." is done! Took $time seconds.";
  return();
}