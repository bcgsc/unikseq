#!/usr/bin/env perl
# RLW
# 11/2020 mod. 11/2021

use strict;

my $version = "v0.2.2 beta";
my ($regsz, $prop, $minnotunique, $minpercentunique) = (100,25,1,90);

if($#ARGV<3){
   print "Usage: $0 $version\n";
   print " < k >\n";
   print " < reference FASTA >\n";
   print " < ingroup FASTA (1 or multi) >\n";
   print " < outgroup FASTA (multi) >\n";
   print " < min. region size (bp) to output (optional, default=$regsz bp) >\n";
   print " < min. average proportion within ingroup (optional, default=$prop %) >\n";
   print " < min. number of non-unique kmer positions allowed in a row (optional, default=$minnotunique) >\n";
   die   " < min. percent unique bases in regions (optional, default=$minpercentunique %) >\n";

}

my $k = $ARGV[0];

my $f1 = $ARGV[1]; #reference
my $f2 = $ARGV[2]; #ingroup
my $f3 = $ARGV[3]; #outgroup

$regsz = $ARGV[4] if($ARGV[4] ne "");
$prop = $ARGV[5] if($ARGV[5] ne "");

$minnotunique = $ARGV[6] if($ARGV[6] ne "");
$minpercentunique = $ARGV[7] if($ARGV[7] ne "");

print "\nRunning: $0 $k $f1 $f2 $f3 $regsz $prop $minnotunique $minpercentunique";

#-----
print "\n\nRecording IDs in $f1 and $f2 to exclude from $f3 ...\n";
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
$fn .= "REGION" . $regsz . "bp-propspcIN" . $prop . "-propunivsOUT" . $minpercentunique;
my $out=$fn . "-uniqueRegions.fa";

&slide($f1,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique);

print "done.\n\nOutput \"unique\" regions of >=$regsz bp are in $out\n";
die "Output unique $k-mers are in $tsv\n";

exit 0;

#--------------------------------
sub slide{

   my ($f,$ex,$in,$k,$incount,$excount,$prop,$out,$tsv,$minnotunique,$minpercentunique) = @_;

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
            &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique);
         }
         $seq = "";
         $prevhead = $head;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
   &printOutput($seq,$prevhead,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique);

   close OUT;
   close TSV;
}

#--------------------------------
sub printOutput{

   my ($seq,$head,$k,$in,$ex,$incount,$excount,$prop,$minnotunique,$minpercentunique) = @_;

   my ($initial,$unique,$notunique,$sum) = (-1,0,0,0);

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

      printf TSV "$pos\t$kmer\toutgroup\t%.4f\n", $ctexf;

      if($ctexf==0){#### kmer is UNIQUE : absent in outgroup, present at sufficient amounts in ingroup!
         $unique++;
         $notunique = 0;### reset notuniquecount
         printf TSV "$pos\t$kmer\tingroup-unique\t%.4f\n", $ctinf;
         $initial = $pos if($initial==-1); ### only track init if was not tracking
         $sum += $ctin;### for average calculation
      }else{      #### kmer not unique!
         $notunique++;
         $sum += $ctin;### for average calculation
         if($notunique > $minnotunique){ ### absent kmer in a row exceed min threshold
            my $stretch = $pos-$initial; ### calculate seq stretch
            my $perunique = $unique / $stretch *100;
            my $avg = $sum / $stretch;#XXX
            my $avgpropspc = $avg / $incount *100;

            #print "$pos..$initial UNIQUE $unique.. STRETCH $stretch.. SEQID $seqid .. AVG $avgpropspc > $prop\n";
            if($stretch>=$regsz && $perunique >= $minpercentunique && $avgpropspc >= $prop){ ### only output longer regions, with a uniqueness prop equal or above user-defined
               my $uniqueseq=substr($seq,$initial,$stretch);
               my $poslast = $pos - 1;
               my $newhead = $head . "region" . $initial . "-" . $poslast . "_size" . $stretch . "_propspcIN";
               printf OUT ">$newhead%.1f" . "_propunivsOUT%.1f" . "\n$uniqueseq\n", ($avgpropspc,$perunique);
            }
            ###reset counters
            $notunique = 0;
            $unique = 0;
            $initial = -1;
            $sum = 0;
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
