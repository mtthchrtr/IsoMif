 #!/usr/bin/perl
use strict;
use warnings;
#use Parallel::ForkManager; #Call the package to use for forks.

#by Matthieu Chartier 10 juin 2013
#Description
#This program runs a batch of isoMif system calls on the jobs list with the isoMif Parameters

#Default values
my $mifPath="";
my $mifParamF="";
my $mifJobsF="";
my $tag="pipeIsoMif";
my $jobsDir="";
my $nbFiles=0;
my $outDir="";
my $cmdMode=0;
my @mifParam=();
my @cases=();
my @casesHT=();
my @cmds=();
my $batch=6;
my $paramString="";

#Read command line
for(my $i=0; $i<=$#ARGV; $i++){
  #print $ARGV[$i]."\n";
  if($ARGV[$i] eq "-e"){ $mifPath=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-p"){ $mifParamF=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-j"){ $mifJobsF=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-c"){ $cmdMode=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-b"){ $batch=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-t"){ $tag=$tag.$ARGV[$i+1]; }
  if($ARGV[$i] eq "-f"){ $jobsDir=$ARGV[$i+1]."/"; }
  if($ARGV[$i] eq "-x"){ $nbFiles=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-o"){ $outDir=$ARGV[$i+1]; }
  if($ARGV[$i] eq "-h"){
    print "##################\nWelcome to pipeIsoMif\n##################\n";
    print "-e         <path to IsoMif program>\n";
    print "-p         <file with all the isoMif parameters>\n";
    print "-j         <file with the list of isoMif jobs>\n";
    print "-c         <cmd mode [nrg (batch qsub on nrg), print (print cmds only), local (run locally) use -b for nb of forks]>\n";
    print "-b         <Nb of cmds per job file / forks if local>\n";
    print "-t         <tag of this job batch>\n";
    print "-f         <dir where to print the job files>\n";
    print "-x         <nb of output files expected in jobsdir>\n";
    print "-o         <dir of the output files>\n";
    print "-h         <print help menu>\n";
    exit;
  }
}

# &storeParams();

&storeCases();

# &recur(0,"",0);

&runCmds();

sub runCmds{
  if($cmdMode eq "nrg"){
    my $filenb=0;
    my $count=0;
    system("rm ".$jobsDir."/*");
    system("rm ".$jobsDir."/ht/*");
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
      system("qsub $file -o /dev/null -e /dev/null");
      # print "qsub $file\n";
    }
    &areJobsDone("nrg",$tag);
  }elsif($cmdMode eq "ht"){
    my $filenb=0;
    my $count=0;
    system("mkdir ".$jobsDir."/ht/") unless(-d $jobsDir."/ht/");
    system("mkdir ".$jobsDir) unless(-d $jobsDir);
    system("rm ".$jobsDir."/*");
    system("rm ".$jobsDir."/ht/*");
    open OUT, ">".$jobsDir.$filenb.".pbs" or die "cant open".$jobsDir.$filenb.".pbs";
    print OUT "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -N ".$tag.$filenb."\n".$mifPath." -pp ".$jobsDir."ht/".$filenb." ".$paramString;
    open OUTHT, ">".$jobsDir."ht/".$filenb or die "cant open".$jobsDir."ht/".$filenb;
    for(my $i=0; $i<@cmds; $i++){
      if($count==$batch){
        $count=0;
        close OUT;
        system("qsub ".$jobsDir.$filenb.".pbs -o /dev/null -e /dev/null");
        $filenb++;
        open OUT, ">".$jobsDir.$filenb.".pbs" or die "cant open".$jobsDir.$filenb.".pbs";
        print OUT "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -N ".$tag.$filenb."\n".$mifPath." -pp ".$jobsDir."ht/".$filenb." ".$paramString;
        open OUTHT, ">".$jobsDir."ht/".$filenb or die "cant open".$jobsDir."ht/".$filenb;
      }
      print OUTHT "$cmds[$i]";
      $count++;
      print OUTHT "\n" unless($count==$batch);
    }
    close OUT;
    close OUTHT;
    system("qsub ".$jobsDir.$filenb.".pbs -o /dev/null -e /dev/null");

    # my @files=glob($jobsDir."/*");
    # foreach my $file (@files){
    #   next unless($file=~/\.pbs$/);
      
    #   # print "qsub $file\n";
    # }

    &areJobsDone("ht",$tag);
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
  }elsif($cmdMode eq "localnrg"){
    for(my $i=0; $i<@cmds; $i++){
      print "$cmds[$i]\n";
      system($cmds[$i]." > /dev/null 2>&1");
    }
  }elsif($cmdMode eq "mammouth"){
    my $filenb=0;
    my $count=0;
    open OUT, ">".$jobsDir.$filenb.".pbs" or die "cant open".$jobsDir.$filenb.".pbs";
    print OUT "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -N ".$tag.$filenb."\n";
    for(my $i=0; $i<@cmds; $i++){
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

    foreach my $file (glob $jobsDir."*"){
      my $call="bqsub_accumulator -q qwork"."@"."mp2 -l walltime=120:00:00 ".$file;
      system($call);
    }
    my $end_call="echo \"echo end\"  | bqsub_accumulator --submit -q qfbb"."@"."mp2 -l walltime=120:00:00";
    system($end_call);
  }
}

sub storeCases{
  if($cmdMode eq "ht"){
    open IN, "<".$mifJobsF or die "Cant open ".$mifJobsF;
    while(my $line=<IN>){
      chomp($line);
      next if($line=~/^$/);
      if($line=~/^#param/){
        $line=~s/^#param\s+//;
        $paramString=$line;
      }else{ push @cmds, $line; }
    }
    close IN;  
  }else{
    open IN, "<".$mifJobsF or die "Cant open ".$mifJobsF;
    while(my $line=<IN>){
      if($line!~/^$/){
        chomp($line);
        push @cases, $mifPath." ".$line;
      }
    }
    close IN;      
  }
}

# sub recur{
#   my $level=$_[2];
#   $level++;
#   for(my $p=$_[0]; $p<@mifParam; $p++){
#     for(my $i=0; $i<@{$mifParam[$p][1]}; $i++){
#       my $cmd=$_[1]." ".$mifParam[$p][0]." ".$mifParam[$p][1][$i]; 
#       if($p==$#mifParam){
#         if($level==@mifParam){
#           foreach my $c (@cases){
#             push @cmds, $c.$cmd;
#           }
#         }
#       }else{
#         &recur($p+1,$cmd,$level);
#       }
#     }
#   }
# }

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
  my $exitLoopC=0;
  my $getout=0;
  my $time=0;
  print "\nWaiting for jobs to terminate..\n";
  if($sys eq "nrg"){
    while(1){
      sleep 15;
      $time+=15;
      my $string;
      if($nbFiles!=0){
        my @nbb=glob $outDir."*";
        print "need $nbFiles, got ".scalar @nbb."\n";
        if(scalar @nbb >= $nbFiles){
          $getout=1;
        }
      }else{
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
        $string="qstat | egrep '".$_[1]."'| egrep ' C ' | wc -l |";
        open COM, $string or die "cant open qstat grep check";
        while(my $line=<COM>){
          print "Completed: $line\n";
          if($line=~/^0$/){
            $exitLoopC=1;
            last;
          }
        }
        close COM;
        $getout=1 if($exitLoopR==1 && $exitLoopQ==1 && $exitLoopC==1); 
      }
      last if($getout==1);
    }
  }elsif($sys eq "ht"){
    while(1){
      sleep 30;
      $time+=30;
      my $got=0;
      my @nbb=glob $outDir."*";
      foreach my $f (@nbb){
        open IN, "<".$f;
        while(my $line=<IN>){
          chomp($line);
          $got++ if($line=~/^REMARK END$/);
        }
        close IN;
      }
      print "need $nbFiles, got ".$got."\n";
      if($got >= $nbFiles){
        $getout=1;
      }
      last if($getout==1);
    }
  }
  print "Job ".$_[0]." is done! Took $time seconds.";
  return();
}