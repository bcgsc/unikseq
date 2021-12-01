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

 < ingroup FASTA (1 or multi) >   tolerated sequences. Used to find regions unique to a % (see option min. average proportion within ingroup)

 < outgroup FASTA (multi) >   outgroup to query kmers against. Note that input reference and ingroup sequences will be automatically excluded from this set.

 < min. region size (bp) to output (optional, default=100 bp) >   minimum "unique" region size to report.

 < min. average proportion within ingroup (optional, default=25 %) > program tracks the numberof qualifying sequences in the ingroup over the sequence stretch, averages and calculates a proportion of the total entries in the ingroup. Sequences are reported only when that proportion is above the minimum set here. 

 < min. number of non-unique kmer positions allowed in a row (optional, default=1) > this controls the tolerance for non-unique kmers in the outgroup when computing a sequence stretch. If set to 0, the sequence stretch is limited to unique kmers only and the % unique bases of the stretch will be 100% (see next parameter). Set to 1, allows some leniency and tolerates a single non-unique kmer in a row. Set to 2, tolerates up to 2 in a row, etc. When set >0, the %unique bases of a sequence stretch will almost never be 100% unique.

 < min. percent unique bases in regions (optional, default=90 %) > controls for sequence uniqueness. 

</pre>


### Plotting "butterfly" plots
-----------

Refer to example.r and replace these lines:

dfa<-read.table("REF_ALFR.fa_IN_ALFR.fa_OUT_iTrackDNA-Database-020821MJA.fa-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("A. fragilis"), " Mt genome"))



