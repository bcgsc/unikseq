# unikseq
identify unique regions in DNA sequences using a kmer approach

## Rene L. Warren, 2020-2021
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

Program to identify unique kmers in a reference, tolerated in an ingroup, not found in outgroup. Useful for designing PCR primer-probe sets with higher specificity 

### Usage
-----------
<pre>

./unikseq.pl v0.2.1 beta
 <k>
 <reference FASTA>
 <ingroup FASTA (1 or multi)>
 <outgroup FASTA (multi)>
 <min. region size (bp) to output (optional, default=100 bp)>
 <min. proportion within ingroup (optional, default=100 %)>

DESCRIPTION OF OPTIONS:
 <k> = kmer length

 <reference FASTA> = reference FASTA (unique FASTA) analysis is done relative to it

 <ingroup FASTA (1 or multi)> = tolerated sequences. Used to find regions unique to a % (see last option, min. proportion)

 <outgroup FASTA (multi)> = outgroup to query kmers against. Note that input reference and ingroup sequences will be automatically excluded from this set.

 <min. region size (bp) to output (optional, default=100 bp)> = minimum "unique" region size to report.

 <min. proportion within ingroup (optional, default=100 %) = minimum % to report unique regions. This should be set to 100 to identify regions common among the ingroup but with no equivalent in the outgroup.

</pre>


### Plotting "butterfly" plots
-----------

Refer to example.r and replace these lines:

dfa<-read.table("REF_ALFR.fa_IN_ALFR.fa_OUT_iTrackDNA-Database-020821MJA.fa-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("A. fragilis"), " Mt genome"))



