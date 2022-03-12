#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
# unikseq

#SYNOPSIS
# Find unique regions in target, tolerated in ingroup but not/less in sequence outgroup

#DOCUMENTATION
#   README.md distributed with this software
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca

#LICENSE
#unikseq Copyright (c) 2020-2022 British Columbia Cancer Agency Branch.  All rights reserved.
#unikseq is released under the GNU General Public License v3

use strict;
use Getopt::Std;

use vars qw($opt_k $opt_r $opt_i $opt_o $opt_s $opt_p $opt_l $opt_u $opt_m);
getopts('k:r:i:o:p:l:u:s:m:');

my $version = "v0-2-8beta";
my ($k, $regsz, $prop, $minnotunique, $minpercentunique,$maxpercentoutgroup) = (25,100,25,1,90,0);

if(! $opt_r || ! $opt_i || ! $opt_o){
   print "Usage: $0 $version\n";
   print "-" x 5, "input files-----\n";
   print " -r reference FASTA (required)\n";
   print " -i ingroup FASTA (required)\n";
   print " -o outgroup FASTA (required)\n";
   print "-" x 5, "kmer uniqueness filters-----\n";
   print " -k length (option, default: -k $k)\n";
   print " -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l $minnotunique)\n";
   print " -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m $maxpercentoutgroup % [original behaviour])\n";
   print "-" x 5, "output filters-----\n";
   print " -s min. reference region [size] (bp) to output (option, default: -s $regsz bp)\n";
   print " -p min. average [proportion] ingroup entries in regions (option, default: -p $prop %)\n";
   die   " -u min. [% unique] kmers in regions (option, default: -u $minpercentunique %)\n";
}

### Fetch options

my $f1 = $opt_r; #reference
my $f2 = $opt_i; #ingroup
my $f3 = $opt_o; #outgroup

$k = $opt_k if($opt_k);
$regsz = $opt_s if($opt_s);
$prop = $opt_p if($opt_p);
$minnotunique = $opt_l if($opt_l);
$minpercentunique = $opt_u if($opt_u);
$maxpercentoutgroup = $opt_m if($opt_m);


###Prepare output
#-----
my $fn = "unikseq_" . $version . "-r_" . $f1 . "-i_" . $f2 . "-o_" . $f3 . "-k" . $k;
my $tsv= $fn . "-uniqueKmers.tsv";

$fn .= "-s" . $regsz . "-p" . $prop . "-l" . $minnotunique . "-u" . $minpercentunique . "-m" . $maxpercentoutgroup;
my $out=$fn . ".fa";
my $log=$fn . ".log";

open(LOG,">$log") || die "Can't write to $log -- fatal.\n";

my $message = "\nRunning: $0 $version\n\t-k $k\n\t-r $f1\n\t-i $f2\n\t-o $f3\n\t-s $regsz\n\t-p $prop\n\t-l $minnotunique\n\t-u $minpercentunique\n\t-m $maxpercentoutgroup\n";

print $message;
print LOG $message;


###Checking files
#-----
if(! -e $f1){
   die "Invalid file: $f1 -- fatal\n";
}

if(! -e $f2){
   die "Invalid file: $f2 -- fatal\n";
}

if(! -e $f3){
   die "Invalid file: $f3 -- fatal\n";
}

#-----
$message = "\nRecording IDs in $f1 and $f2 to exclude from $f3 ...\n";

print $message;
print LOG $message;

my $rec;
open(IN,$f1) || die "Can't read $f1 -- fatal.\n";
while(<IN>){
   chomp;
   s/\r\n/\n/g;### DOS to UNIX
   $rec->{$1}=1 if(/^\>(\S+)/);
}
close IN;
#-----
open(IN,$f2) || die "Can't read $f2 -- fatal.\n";
while(<IN>){
   chomp;
   s/\r\n/\n/g;### DOS to UNIX
   $rec->{$1}=1 if(/^\>(\S+)/);
}
close IN;
#-----

$message = "done.\nReading outgroup $f3, excluding records in $f1 and $f2 ...\n";
print $message;
print LOG $message;

my ($ex,$excount) = &readFasta($f3,$k,$rec);##exclude outgroup
my $rec; ### re-initialize record

$message = "done.\nReading ingroup $f2 ...\n";
print $message;
print LOG $message;

my ($in,$incount) = &readFasta($f2,$k,$rec);##include ingroup

$message = "done.\nBeginning kmer analysis (k$k), sliding base by base on $f1 ...\n";
print $message;
print LOG $message;

&slide($f1,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique,$maxpercentoutgroup);

$message = "done.\n";
$message .= "-" x 30, "\n";
$message .= "\nOutput unique reference sequence regions >= $regsz bp in:\n$out\n";
$message .= "\nOutput unique $k-mers (for butterfly plot):\n$tsv\n\n";

print $message;
print LOG $message;

close LOG;

exit;

#--------------------------------
sub slide{

   my ($f,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique,$maxpercentoutgroup) = @_;

   my ($head,$prevhead,$seq) = ("","","");

   open(OUT,">$out") || die "Can't write $out -- fatal.\n";
   open(TSV,">$tsv") || die "Can't write $tsv -- fatal.\n";

   print TSV "position\tkmer\tcondition\tvalue\n";

   open(IN,$f) || die "Can't read $f -- fatal.\n";
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){
         $head = $1;
         if($prevhead ne $head && $prevhead ne "" && $seq ne ""){
            &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup);
         }
         $seq = "";
         $prevhead = $head;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup);

   close OUT;
   close TSV;
}

#--------------------------------
sub printOutput{

   my ($seq,$head,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup) = @_;

   my ($initial,$unique,$notunique,$sum,$sumout) = (-1,0,0,0,0);

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $kmer = uc($kmer);
      #print "$kmer ......\n";

      my $listex = $ex->{$kmer};
      my $ctex = keys(%$listex);
      my $listin = $in->{$kmer};
      my $ctin = keys(%$listin);

      my $ctexf = $ctex/$excount;##as a fraction
      my $ctinf = $ctin/$incount;##as a fraction
      $ctinf = -1 * $ctinf;

      #print "$ctin .. $incount $ctinf %%\n";

      if($ctexf==0 || ($ctexf*100) < $maxpercentoutgroup){#### kmer is UNIQUE : absent in outgroup, present at sufficient amounts in ingroup!

         $unique++;
         $notunique = 0;### reset notuniquecount
         printf TSV "$pos\t$kmer\tingroup-unique\t%.4f\n", $ctinf;
         $initial = $pos if($initial==-1); ### only track init if was not tracking
         $sum += $ctin;### for average calculation
         $sumout += $ctex;

      }else{      #### kmer not unique!

         printf TSV "$pos\t$kmer\toutgroup\t%.4f\n", $ctexf; 

         if($initial > -1){###do not update trackers unless the region started with unique seqs. 
            $notunique++;
            $sum += $ctin;### for average calculation

            if($notunique > $minnotunique){ ### absent kmer in a row exceed min allowed threshold
               $sum -= $ctin;
               my $stretch = $pos-$initial; ### calculate seq stretch
               my $perunique = $unique / $stretch *100;
               my $avg = $sum / $stretch;#XXX
               my $avgout = $sumout / $stretch;
               my $avgpropspc = $avg / $incount *100;

               #print "$pos..$initial UNIQUE $unique.. STRETCH $stretch.. SEQID $seqid .. AVG $avgpropspc > $prop\n";
               if($stretch>=$regsz && $perunique >= $minpercentunique && $avgpropspc >= $prop){ ### only output longer regions, with a uniqueness prop equal or above user-defined
                  my $uniqueseq=substr($seq,$initial,$stretch);
                  my $poslast = $pos - 1;
                  my $newhead = $head . "region" . $initial . "-" . $poslast . "_size" . $stretch . "_propspcIN";
                  printf OUT ">$newhead%.1f" . "_propunivsOUT%.1f_avgOUTentries%.1f" . "\n$uniqueseq\n", ($avgpropspc,$perunique,$avgout);
               }
               ###reset counters
               $notunique = 0;
               $unique = 0;
               $initial = -1;
               $sum = 0;
               $sumout = 0;
            }
         }
      }
   }
}

#--------------------------------
sub readFasta{

   my ($f,$k,$rec) = @_;
   my $h;
   my ($head,$prevhead,$seq) = ("","","");
   my ($flag,$prevflag,$count) = (0,0,0);
   
   open(IN,$f) || die "Can't read $f -- fatal.\n";
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){

         $head = $1;

         if(defined $rec->{$head}){
            $flag=1;
         }else{
            $flag=0;
            $count++;
         }

         if ($head ne $prevhead && $seq ne '' && $prevhead ne '' && $prevflag==0){
            $h = &kmerize($seq,$k,$prevhead,$h);
         }
         $seq = '';
         $prevhead = $head;
         $prevflag = $flag;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   if ($prevflag==0){
      $h = &kmerize($seq,$k,$prevhead,$h);
   }

   close IN;
   return $h,$count;
}

#-----------------------
sub kmerize{

   my ($seq,$k,$head,$h) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $kmer = uc($kmer);
      my $rckmer = reverseComplement($kmer);
      #print "$kmer .. $rckmer\n";
      $h->{$kmer}{$head}++;
      $h->{$rckmer}{$head}++;
   }
   return $h;
}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/;
   return (reverse());
}
