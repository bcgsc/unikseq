#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
#writeBloom

#SYNOPSIS

#DOCUMENTATION
#   This code was originally writtent for PERL LINKS. Please see this readme for details:
#   https://github.com/bcgsc/LINKS#testing-the-bloom-filters
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use LINKS, unikseq, the LINKS/unikseq code or ideas, please cite our work

#LICENSE
#   LINKS Copyright (c) 2014-2023 Canada's Michael Smith Genome Science Centre.  All rights reserved.
#   Unikseq Copyright (c) 2020-2023 Canada's Michael Smith Genome Science Centre.  All rights reserved.
#   Unikseq is released under the GNU General Public License v3
#   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.
#   For commercial licensing options, please contact Patrick Rebstein prebstein@bccancer.bc.ca

use strict;
use POSIX;
use FindBin;
use lib "$FindBin::RealBin/./libexec";
use BloomFilter;
use Getopt::Std;
use Net::SMTP;
use Time::HiRes;
use vars qw($opt_f $opt_k $opt_p);
getopts('f:k:p:');
my ($k,$fpr)=(25,0.001);

#-------------------------------------------------

if(! $opt_f ){
   print "Usage: $0\n";
   print "-f  sequences to scaffold (Multi-FASTA format, required)\n"; 
   print "-k  k-mer value (default -k $k, optional)\n";
   die "-p  Bloom filter false positive rate (default -p $fpr, optional - increase to prevent memory allocation errors)\n";
}

my $assemblyfile = $opt_f;
$k = $opt_k if($opt_k);
$fpr = $opt_p if($opt_p); 

print "\nRunning:$0 -f $assemblyfile -k $k -p $fpr\n\n";

my $bfout = $assemblyfile . "_k" . $k . "_p" . $fpr . "_rolling.bloom";

if(! -e $assemblyfile){
   my $file_message = "\nInvalid file: $assemblyfile -- fatal\n";
   print $file_message;
   exit;
}else{
   my $file_message = "Checking sequence target file $assemblyfile...ok\n";
   print $file_message;
}


my $date = `date`;
chomp($date);
my $bloom;

eval{


print "$date:Estimating number of elements from file size\n";

my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = stat($assemblyfile);
my $bfelements = int($size);

my $date = `date`;
chomp($date);

my $m = ceil((-1 * $bfelements * log($fpr)) / (log(2) * log(2))); #the number of bits
my $rem = 64 - ($m % 64);
$m = $m + $rem;

my $hashfct = floor(($m / $bfelements) * log(2));# the number of hash functions
print "*****\nBloom filter specs\nelements=$bfelements\nFPR=$fpr\nsize (bits)=$m\nhash functions=$hashfct\n*****\n";
$bloom = new BloomFilter::BloomFilter($m, $hashfct, $k);

print "$date:Shredding supplied sequence file (-f $assemblyfile) into $k-mers..\n";
$bloom = &contigsToBloom($assemblyfile,$hashfct,$k,$bloom);

my $date = `date`;
chomp($date);


my $writing_tigbloom_message = "\n$date:Writing Bloom filter to disk ($bfout)\n";
print $writing_tigbloom_message;

$bloom->storeFilter($bfout);

};

my $date = `date`;
chomp($date);

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
}else{
   my $success = "\n$date:$0 executed normally\n";
   print $success;
}

# manually free the data structures 
#$bloom->DESTROY;

exit 1;


#----------------
sub contigsToBloom{
   my ($file,$hashfct,$k,$bloom) = @_;

   my $prevhead = "";
   my $seq = "";
   my $cttig=0;
   open(IN,$file) || die "Error reading $file -- fatal.\n";

   print "Contigs processed k=$k:\n";
   ###
   while(<IN>){
      chomp;
      if(/^\>(\S+)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;
            $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k);
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   $|++;
   $bloom = &kmerizeContigBloom_newloop(uc($seq),$bloom,$hashfct,$k);
   close IN;
   
   return $bloom;
}

#----------------
sub kmerizeContigBloom{
   my ($seq,$bloom,$hashfuct, $k) = @_;

   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      $bloom->insert($kmer);
   }
   return $bloom;
}

#-----------------------
sub kmerizeContigBloom_new{
    my ($seq,$bloom,$hashfct,$k) = @_;

    my $itr = new BloomFilter::RollingHashIterator($seq, $hashfct, $k);
    my $next = $itr->getNext();
    while($next) {
        $bloom->insert($next);
        $next = $itr->getNext();
    }
    return $bloom;
}

#-----------------------
sub kmerizeContigBloom_newloop{
    my ($seq,$bloom,$hashfct,$k) = @_;
    
    BloomFilter::insertSeq($bloom, $seq, $hashfct, $k);
    return $bloom
}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/; 
   return (reverse());
}
            
