# unikseq
identify unique regions in DNA sequences using a kmer approach

## Rene L. Warren, 2020-2021
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

Program to identify unique kmers in a reference, tolerated in an ingroup, not found in outgroup. Useful for designing PCR primer-probe sets with higher specificity 

### Running unikseq
-----------

<pre>
Usage: ./unikseq.pl v0.2.2 beta
 < k >
 < reference FASTA >
 < ingroup FASTA (1 or multi) >
 < outgroup FASTA (multi) >
 < min. region size (bp) to output (optional, default=100 bp) >
 < min. average proportion within ingroup (optional, default=25 %) >
 < min. number of non-unique kmer positions allowed in a row (optional, default=1) >
 < min. percent unique bases in regions (optional, default=90 %) >
</pre>

Notes:
<pre>

 < k >   kmer length

 < reference FASTA >   reference FASTA (unique FASTA) analysis is done relative to it

 < ingroup FASTA (1 or multi) >   tolerated sequences. Used to find regions unique to a % (see last option, min. proportion)

 < outgroup FASTA (multi) >   outgroup to query kmers against. Note that input reference and ingroup sequences will be automatically excluded from this set.

 < min. region size (bp) to output (optional, default=100 bp) >   minimum "unique" region size to report.

 < min. proportion within ingroup (optional, default=100 %) >   minimum % to report unique regions. This should be set to 100 to identify regions common among the ingroup but with no equivalent in the outgroup.

</pre>


### Plotting "butterfly" plots
-----------

Refer to example.r and replace these lines:

dfa<-read.table("REF_ALFR.fa_IN_ALFR.fa_OUT_iTrackDNA-Database-020821MJA.fa-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("A. fragilis"), " Mt genome"))



