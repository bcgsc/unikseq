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
Usage: ./unikseq.pl [v0.2.4 beta]
 -r reference FASTA (required)
 -i ingroup FASTA (required)
 -o outgroup FASTA (required)
 -k length (option, default: -k 25)
 -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
 -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
 -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
 -u min. [% unique] kmers in regions (option, default: -u 90 %)
</pre>

Notes:
<pre>

 -k length (option, default: -k 25)
  kmer length

 -r reference FASTA (required)
  reference FASTA (unique FASTA) analysis is done relative to it

 -i ingroup FASTA (required)
  tolerated sequences. Used to find regions unique to a % (see -p option)

 -o outgroup FASTA (required)
  outgroup to query kmers against. Note that input reference and ingroup sequences will be automatically excluded from this set.

 -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
  minimum "unique" reference (target) region size to report.

 -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
  program tracks the number of qualifying sequences in the ingroup over the sequence stretch, averages and calculates a proportion of the total entries in the ingroup. Sequences are reported only when that proportion is above the minimum set here. 

 -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
  this leniency factor controls the tolerance for non-unique kmers in the outgroup when computing a sequence stretch. If set to 0, the sequence stretch is limited to unique kmers only and the % unique bases of the stretch will be 100% (see next parameter). Set to 1, allows some leniency and tolerates a single non-unique kmer in a row. Set to 2, tolerates up to 2 in a row, etc. When set >0, the %unique bases of a sequence stretch will almost never be 100% unique and so -u MUST be lowered below 100%.

 -u min. [% unique] kmers in regions (option, default: -u 90 %)
  controls for sequence uniqueness in the reference output regions. 

</pre>


### Plotting "butterfly" plots
-----------

Refer to example.r and replace these lines:

dfa<-read.table("REF_ALFR.fa_IN_ALFR.fa_OUT_iTrackDNA-Database-020821MJA.fa-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("A. fragilis"), " Mt genome"))



