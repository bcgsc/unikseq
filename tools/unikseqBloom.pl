#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
# unikseqBloom - a simple prototype support for large genomes, with Bloom filter functionality

#SYNOPSIS
# Find unique regions in target, tolerated in ingroup but not/less in sequence outgroup

#DOCUMENTATION
#   README.md distributed with this software
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca

#LICENSE
#   Unikseq Copyright (c) 2020-2023 British Columbia Cancer Agency Branch.  All rights reserved.
#   Unikseq is released under the GNU General Public License v3
#   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.
#   For commercial licensing options, please contact Patrick Rebstein prebstein@bccancer.bc.ca

use strict;
use Getopt::Std;
###
use POSIX;
use FindBin;
use lib "$FindBin::RealBin/./libexec";
use BloomFilter;
###
use vars qw($opt_k $opt_r $opt_i $opt_o $opt_s $opt_p $opt_l $opt_u $opt_m $opt_c $opt_t $opt_v);
getopts('k:r:i:o:p:l:u:s:m:c:t:v:');

my $version = "v1.0.0";
my ($k, $regsz, $prop, $minnotunique, $minpercentunique,$maxpercentoutgroup,$cflag,$tsvflag) = (25,100,0,0,90,0,0,0);

if(! $opt_r || ! $opt_i || ! $opt_o){
   print "Usage: $0 $version\n";
   print "-" x 5, "input files-----\n";
   print " -r reference FASTA (required)\n";
   print " -i ingroup Bloom filter (required)\n";
   print " -o outgroup Bloom filter (required)\n";
   print "-" x 5, "kmer uniqueness filters-----\n";
   print " -k length (option, default: -k $k)\n";
   print " -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l $minnotunique)\n";
   print " -m max. [% entries] in outgroup tolerated to have a reference kmer (option, default: -m $maxpercentoutgroup % [original behaviour])\n";
   print "-" x 5, "output filters-----\n";
   print " -t print only first t bases in tsv output (option, default: -t [k])\n";
   print " -c output conserved FASTA regions between reference and ingroup entries (option, -c 1==yes -c $cflag==no [default, original unikseq behaviour])\n";
   print " -s min. reference FASTA region [size] (bp) to output (option, default: -s $regsz bp)\n";
   print " -p min. [-c 0:region average /-c 1: per position] proportion of ingroup entries (option, default: -p $prop %)\n";
   print " -u min. [% unique] kmers in regions (option, default: -u $minpercentunique %)\n";
   die " -v output tsv files (option, -v 1==yes -v 0==no [default])\n";
}

### Fetch options

my $f1 = $opt_r; #reference
my $f2 = $opt_i; #ingroup
my $f3 = $opt_o; #outgroup

$k = $opt_k if($opt_k);
my $tchar = $k;
$regsz = $opt_s if($opt_s);
$prop = $opt_p if($opt_p);
$cflag = $opt_c if($opt_c);
$minnotunique = $opt_l if($opt_l);
$minpercentunique = $opt_u if($opt_u);
$maxpercentoutgroup = $opt_m if($opt_m);
$tchar = $opt_t if($opt_t);
$tsvflag = $opt_v if($opt_v);
my $MAXSTRINGLEN = 5000000;


###Prepare output
#-----
my $fn = "unikseq_" . $version . "-r_" . $f1 . "-i_" . $f2 . "-o_" . $f3 . "-k" . $k;
my $tsv= $fn . "-uniqueKmers.tsv";
my $tsvcons = $fn . "-conservedKmers.tsv";

$fn .= "-c" . $cflag . "-s" . $regsz . "-p" . $prop . "-l" . $minnotunique . "-u" . $minpercentunique . "-m" . $maxpercentoutgroup;

my $outunique=$fn . "-unique.fa";
my $outcons=$fn . "-conserved.fa";

my $log=$fn . ".log";

open(LOG,">$log") || die "Can't write to $log -- fatal.\n";

my $message = "\nRunning: $0 $version\n\t-k $k\n\t-r $f1\n\t-i $f2\n\t-o $f3\n\t-t $tchar\n\t-c $cflag\n\t-s $regsz\n\t-p $prop\n\t-l $minnotunique\n\t-u $minpercentunique\n\t-m $maxpercentoutgroup\n\t-v $tsvflag\n";

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

if($tchar < 1 || $tchar > $k){
   die "Option -t must be a valid number between 1 and $k -- fatal\n";
}

#-----

$message = "done.\nReading outgroup $f3 ...\n";
print $message;
print LOG $message;

my $excount = 1;
#my ($ex,$excount) = &readFasta($f3,$k);##exclude outgroup
my $ex = new BloomFilter::BloomFilter($f3);

$message = "done.\nReading ingroup $f2 ...\n";
print $message;
print LOG $message;

my $incount = 1;
#my ($in,$incount) = &readFasta($f2,$k);##include ingroup
my $in = new BloomFilter::BloomFilter($f2);

$message = "done.\nBeginning kmer analysis (k$k), sliding base by base on $f1 ...\n";
print $message;
print LOG $message;

my $conseq;

if($cflag){
   $conseq = &slideConserved($f1,$ex,$in,$k,$incount,$excount,$prop,$outcons,$tsvcons,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag);
   $message = "\ndone.\n";
   $message .= "-" x 30, "\n";
   $message .= "\nOutput conserved reference sequence regions >= $regsz bp (vs. ingroup FASTA) in:\n$outcons\n";

   $message .= "\nOutput conserved $k-mers within ingroup FASTA (for your reference):\n$tsvcons\n\n";
   print $message;
   print LOG $message;
}

&slide($cflag,$conseq,$f1,$ex,$in,$k,$incount,$excount,$prop,$outunique,$tsv,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag);

$message = "\ndone.\n";
$message .= "-" x 30, "\n";
$message .= "\nOutput unique reference sequence regions >= $regsz bp in:\n$outunique\n";
$message .= "\nOutput unique $k-mers (for butterfly plot):\n$tsv\ndone.\n";
print $message;
print LOG $message;

close LOG;

exit;

#--------------------------------
sub slideConserved{

   my ($f,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag) = @_;

   my ($head,$prevhead,$seq,$preventry) = ("","","","");
   my $conseq;

   open(OUT,">$out") || die "Can't write $out -- fatal.\n";

   if($tsvflag){
      open(TSV,">$tsv") || die "Can't write $tsv -- fatal.\n";
      print TSV "header\tposition\t$tchar-mer\tcondition\tnum_entries\tproportion\n";
   }

   my $ctseq=0;
   open(IN,$f) || die "Can't read $f -- fatal.\n";
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){

         $ctseq++;
         print "\r$ctseq";
         $|++;

         my $entry=$1;
         $head = $_;

         if($prevhead ne $head && $prevhead ne "" && $seq ne ""){
            $conseq = &printConserved($conseq,$seq,$preventry,$k,$in,$incount,$prop,$tchar,$tsvflag);
         }
         $seq = "";
         $prevhead = $head;
         $preventry = $entry;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   $conseq = &printConserved($conseq,$seq,$preventry,$k,$in,$incount,$prop,$tchar,$tsvflag);

   close OUT;
   close TSV if($tsvflag);

   return $conseq;
}

#--------------------------------
sub slide{

   my ($cflag,$conseq,$f,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag) = @_;

   my ($head,$prevhead,$seq,$preventry) = ("","","","");

   open(OUT,">$out") || die "Can't write $out -- fatal.\n";

   if($tsvflag){
      open(TSV,">$tsv") || die "Can't write $tsv -- fatal.\n";
      print TSV "header\tposition\t$tchar-mer\tcondition\tvalue\n";
   }

   open(IN,$f) || die "Can't read $f -- fatal.\n";
   my $ctseq = 0;
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){

         $ctseq++;
         print "\r$ctseq";
         $|++;

         my $entry = $1;
         $head = $_;
         if($prevhead ne $head && $prevhead ne "" && $seq ne ""){
            &printOutput($cflag,$conseq,$seq,$preventry,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag);
         }
         $seq = "";
         $prevhead = $head;
         $preventry = $entry;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   &printOutput($cflag,$conseq,$seq,$preventry,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag);

   close OUT;
   close TSV if($tsvflag);
}

#--------------------------------
sub printConserved{

   my ($conseq,$masterseq,$head,$k,$in,$incount,$prop,$tchar,$tsvflag) = @_;

   ### sequence/sub-specific variables
   my $newmaster = "";
   my $last2k = "";
   my $adjustedpos = 0;
   my $mastercoord = 0;
   my $chunknum = 0;

   my $consequ="";

   print "\nprocessing chunks...\n";

   for(my $mpos=0; $mpos<= length($masterseq); $mpos+=$MAXSTRINGLEN){
      $chunknum++;
      my $seq = substr($masterseq,$mastercoord,$MAXSTRINGLEN);###This is needed for better memory management on large strings
      my $tmpconsequ = lc($seq);
      my $len = length($seq);
      print "Adjusted position: $adjustedpos Tile-mpos:$mpos, mastercoord: $mastercoord step: $MAXSTRINGLEN (length of seq = $len)\n" if($tsvflag);

      my ($initial,$sum,$buffer) = (-1,0,0);

      POSITION:
      for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
         my $kmer = substr($seq,$pos,$k);
         my $tcharmer = substr($kmer,0,$tchar);
         $kmer =~ tr/U/T/; ### handles RNA U>>>T

         #my $listin = $in->{$kmer};
         #my $ctin = keys(%$listin);

         my ($ctexf,$ctinf) = (0,0);
         $ctinf = 1 if($in->contains($kmer));    # $ctin/$incount if($incount);##as a fraction
         if($tsvflag){
            printf TSV "$head\t$pos\t$tcharmer\tingroup\t$ctinf\t%.4f\n", $ctinf if($tsvflag);
         }
         if((abs($ctinf)*100) >= $prop){#### kmer is conserved at set proportion in ingroup

            $initial = $pos if($initial==-1); ### only track init if was not tracking
            $sum += $ctinf;### for average calculation

            if($pos==(length($seq)-$k)){###end of sequence
               $buffer = 1;
               my ($adjstart,$adjend)=($initial+$adjustedpos,$pos+$adjustedpos);
               ($initial,$sum,$tmpconsequ) = &outputFASTA($initial,$ctinf,$sum,$prop,$regsz,$head,$seq,$pos,$tmpconsequ,$buffer,$adjstart,$adjend);
            }         
         }else{#### kmer below set threshold
               my ($adjstart,$adjend)=($initial+$adjustedpos,$pos+$adjustedpos);
               ($initial,$sum,$tmpconsequ) = &outputFASTA($initial,$ctinf,$sum,$prop,$regsz,$head,$seq,$pos,$tmpconsequ,$buffer,$adjstart,$adjend);
         }
      }
   
      my $seqminuslastk = substr($tmpconsequ, 0, (length($tmpconsequ) - (2*$k)));###The last kmer can never be resolved (in fact last <2*k can's be resolved)
      $last2k = substr($tmpconsequ, (length($tmpconsequ) - (2*$k)), (2*$k));
      $newmaster .= $seqminuslastk;
      $adjustedpos = length($newmaster);
      $mastercoord = ($mpos + $MAXSTRINGLEN) - ((2*$k)*$chunknum);

   }###master seq loop

   $newmaster .= $last2k;

   $conseq->{$head} = $newmaster;

   return $conseq;
}

#--------------------------------
sub outputFASTA{

   my ($initial,$ctin,$sum,$prop,$regsz,$head,$seq,$pos,$consequ,$buffer,$adjstart,$adjend) = @_;

   if($initial > -1){###do not update trackers unless the region started with conserved seqs.
      my $stretch = $pos-$initial+$buffer; ### calculate seq stretch
      my $avg = $sum / $stretch;
      my $avgpropspc = $avg / $incount *100;

      my $conservedseq=substr($seq,$initial,$stretch);
      substr($consequ, $initial, $stretch, uc($conservedseq)); ### !!! NOT EFFICIENT IN PERL, AND MAY NOT WORK ON LARGE STRINGS -- see ntEdit for solution

      if($stretch>=$regsz && $avgpropspc >= $prop){ ### only output longer regions, with a ingroup prop equal or above user-defined
         my $poslast = $adjend - 1 + $buffer;
         my $newhead = $head . "region" . $adjstart . "-" . $poslast . "_size" . $stretch . "_propspcIN";
         printf OUT ">$newhead%.1f" . "\n$conservedseq\n", ($avgpropspc);
      }
      ###reset counters
      $initial = -1;
      $sum = 0;
   }
   return $initial,$sum,$consequ;
}

#--------------------------------
sub printOutput{

   my ($cflag,$conseq,$seq,$head,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique,$maxpercentoutgroup,$tchar,$tsvflag) = @_;

   my ($initial,$unique,$notunique,$sum,$sumout,$buffer) = (-1,0,0,0,0,0);

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      my $tcharmer = substr($kmer,0,$tchar);
      $kmer =~ tr/U/T/; ### handles RNA U>>>T

      ###my $listex = $ex->{$kmer};
      ###my $ctex = keys(%$listex);
      #my $listin = $in->{$kmer};
      #my $ctin = keys(%$listin);

      my ($ctexf,$ctinf) = (0,0);

      $ctexf = 1 if($ex->contains($kmer)); #$ctex/$excount if($excount);##as a fraction
      $ctinf = 1 if($in->contains($kmer)); #$ctin/$incount if($incount);##as a fraction
      $ctinf = -1 * $ctinf;

      if(($pos < (length($seq)-$k)) && ($ctexf==0 || ($ctexf*100) < $maxpercentoutgroup)){#### kmer is UNIQUE : absent in outgroup, present at sufficient amounts in ingroup!

         $unique++;
         $notunique = 0;### reset notuniquecount
         if($tsvflag){
            printf TSV "$head\t$pos\t$tcharmer\tingroup-unique\t%.4f\n", $ctinf if($tsvflag);
         }
         $initial = $pos if($initial==-1); ### only track init if was not tracking
         $sum += abs($ctinf);### for average calculation
         $sumout += $ctexf;

      }else{      #### kmer not unique!

         if($ctexf){
            printf TSV "$head\t$pos\t$tcharmer\toutgroup\t%.4f\n", $ctexf if($tsvflag);
         }elsif($pos == (length($seq)-$k) && $ctinf<=0){###end of sequence exception
            $unique++;
            $notunique = 0;### reset notuniquecount
            printf TSV "$head\t$pos\t$tcharmer\tingroup-unique\t%.4f\n", $ctinf if($tsvflag);
            $sum += abs($ctinf);### for average calculation
            $sumout += $ctexf;
            $buffer = 1;
         }

         if($initial > -1){###do not update trackers unless the region started with unique seqs. 
            $notunique++;
            $sum += $incount;### for average calculation

            if($notunique > $minnotunique){ ### absent kmer in a row exceed min allowed threshold
               $sum -= $incount;
               my $stretch = $pos-$initial+$buffer; ### calculate seq stretch
               my $perunique = $unique / $stretch *100;
               my $avg = $sum / $stretch;
               my $avgout = $sumout / $stretch;
               my $avgpropspc = $avg / $incount *100;

               if($cflag){
                  if($stretch>=$regsz && $perunique >= $minpercentunique){
                     my $uniqueseq=substr($conseq->{$head},$initial,$stretch);
                     my $poslast = $pos - 1 + $buffer;
                     my $newhead = $head . "region" . $initial . "-" . $poslast . "_size" . $stretch . "_propspcIN";
                     printf OUT ">$newhead%.1f" . "_propunivsOUT%.1f_avgOUTentries%.1f" . "\n$uniqueseq\n", ($avgpropspc,$perunique,$avgout);
                  }
               }else{#original behaviour

#print "$stretch>=$regsz && $perunique >= $minpercentunique && $avgpropspc >= $prop ****\n";

                  if($stretch>=$regsz && $perunique >= $minpercentunique && $avgpropspc >= $prop){ ### only output longer regions, with a uniqueness prop equal or above user-defined
                     my $uniqueseq=substr($seq,$initial,$stretch);
                     my $poslast = $pos - 1 + $buffer;
                     my $newhead = $head . "region" . $initial . "-" . $poslast . "_size" . $stretch . "_propspcIN";
                     printf OUT ">$newhead%.1f" . "_propunivsOUT%.1f_avgOUTentries%.1f" . "\n$uniqueseq\n", ($avgpropspc,$perunique,$avgout);
                  }
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

   my ($f,$k) = @_;

   my $h;
   my $rec;
   my ($head,$prevhead,$seq,$preventry) = ("","","","");
   
   open(IN,$f) || die "Can't read $f -- fatal.\n";
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){### first string with non-space characters is unique ID for entry/entries
         my $entry=$1;
         $head = $_;
         $rec->{$entry} = 1 ;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $h = &kmerize($seq,$k,$preventry,$h);
         }

         $seq = '';
         $prevhead = $_;
         $preventry = $entry;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }

   $h = &kmerize($seq,$k,$preventry,$h);

   close IN;

   my $count = keys(%$rec);###added 2022.10.08
   print "***$count unique entries***\n";

   return $h,$count;
}

#-----------------------
sub kmerize{

   my ($seq,$k,$head,$h) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $kmer =~ tr/U/T/; ### handles RNA U>>>T
      my $rckmer = reverseComplement($kmer);
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
