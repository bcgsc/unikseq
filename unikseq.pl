#!/usr/bin/env perl
# RLW
# 11/2020 mod. 11/2021

use strict;

my $version = "v0.2.1 beta";
my $regsz = 100;
my $prop = 100;

if($#ARGV<3){
   print "Usage: $0 $version\n";
   print " <k>\n";
   print " <reference FASTA>\n";
   print " <ingroup FASTA (1 or multi)>\n";
   print " <outgroup FASTA (multi)>\n";
   print " <min. region size (bp) to output (optional, default=$regsz bp)>\n";
   die " <min. proportion within ingroup (optional, default=$prop %)>\n";
}

my $k = $ARGV[0];

my $f1 = $ARGV[1]; #reference
my $f2 = $ARGV[2]; #ingroup
my $f3 = $ARGV[3]; #outgroup

$regsz = $ARGV[4] if($ARGV[4] ne "");
$prop = $ARGV[5] if($ARGV[5] ne "");

print "Running: $0 $k $f1 $f2 $f3 $regsz $prop";

#-----
print "\nRecording IDs in $f1 and $f2 to exclude from $f3 ...\n";
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

print "done.\nReading outgroup $f3, excluding records in $f1 and $f2 ...\n";
my ($ex,$excount) = &readFasta($f3,$k,$rec);##exclude outgroup
my $rec; ### re-initialize record
print "done.\nReading ingroup $f2...\n";
my ($in,$incount) = &readFasta($f2,$k,$rec);##include ingroup

print "done.\nBeginning kmer analysis (k$k), sliding base by base on $f1 ...\n";

my $fn = "REF_" . $f1 . "_IN_" . $f2 . "_OUT_" . $f3;
my $tsv=$fn . "-uniqueKmers.tsv";
$fn .= "REGION" . $regsz . "bp-PROP" . $prop;
my $out=$fn . "-uniqueRegions.fa";

&slide($f1,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv);

print "done.\n\nOutput \"unique\" regions of >=$regsz bp are in $out\n";
die "Output unique $k-mers are in $tsv\n";

exit 0;

#--------------------------------
sub slide{

   my ($f,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv) = @_;

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
            &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop);
         }
         $seq = "";
         $prevhead = $head;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop);

   close OUT;
   close TSV;
}

#--------------------------------
sub printOutput{

   my ($seq,$head,$k,$in,$ex,$incount,$excount,$prop) = @_;

   my $initial = 0;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $kmer = uc($kmer);
      #print "$kmer ......\n";

      my $listex = $ex->{$kmer};
      my $ctex = keys(%$listex);
      my $listin = $in->{$kmer};
      my $ctin = keys(%$listin);

      $ctex = $ctex/$excount;##as a fraction
      $ctin = $ctin/$incount;##as a fraction
      $ctin = -1 * $ctin;

      printf TSV "$pos\t$kmer\toutgroup\t%.4f\n", $ctex;

      if(! $ctex){#### kmer is UNIQUE!
         printf TSV "$pos\t$kmer\tingroup-unique\t%.4f\n", $ctin;
         $initial = $pos if($initial==0); ### only track init if was not tracking
      }else{      #### kmer not unique!
         my $stretch=$pos-$initial+1; ### calculate seq stretch
         if($initial && $stretch>=$regsz && (abs($ctin)*100) >= $prop ){ ### only output longer regions, with a uniqueness prop equal or above user-defined
            my $uniqueseq=substr($seq,$initial,$stretch);
            my $actprop = abs($ctin)*100;
            my $newhead = $head . "region" . $initial . "-" . $pos . "_size" . $stretch . "_propIngroupSpecies";
            printf OUT ">$newhead%.1f\n$uniqueseq\n",$actprop;
         }
         $initial=0;
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
