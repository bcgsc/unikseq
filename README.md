# unikseq
identify unique regions in DNA sequences using a kmer approach

## Rene L. Warren, 2020-2022
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

Program to identify unique kmers in a reference, tolerated in an ingroup, not found in outgroup. Useful for designing PCR primer-probe sets with higher specificity 

### Running unikseq
-----------

<pre>
Usage: ./unikseq.pl v0-2-8beta 
 -r reference FASTA (required)
 -i ingroup FASTA (required)
 -o outgroup FASTA (required)
 -k length (option, default: -k 25)
 -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
 -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
 -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
 -u min. [% unique] kmers in regions (option, default: -u 90 %)
 -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m 0 % [original behaviour])
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

 -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m 0 % [original behaviour])
  controls the kmer "uniqueness" in the outgroup, by tolerating a certain fraction of sequence entries in the outgroup having the reference kmer. This option could be useful when there's high similarity between the reference, ingroup AND outgroup sequences and more fine-grain adjustments are needed. Not specifying this option (-m 0) is the original unikseq behaviour.

 Example command:
 ./unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90

</pre>


### Outputs (2)

1) TSV file (.tsv) 

   Tab-Separated Variable file. Reports all reference sequence kmers in 4 columns:
   <pre>
   position	kmer	condition	value
   [coordinates][sequence][in/out group][proportion (rate) in each in/out group]
   By default, every instance of a reference kmer is reported when found in the outgroup.
   When it is not found, the ingroup-unique will be reported (if found).
   If a reference kmer is found in outgroup sequences, the ingroup-unique WILL NOT report any
   values.

   e.g.
   0	GCTAGTGTAGCTTAATGTAAAGTAT	ingroup-unique	-0.1270
   ...
   207	ACCTTGCTAAGCCACACCCCCAAGG	ingroup-unique	-0.9683
   208	CCTTGCTAAGCCACACCCCCAAGGG	outgroup	0.0115
   209	CTTGCTAAGCCACACCCCCAAGGGA	outgroup	0.0115
   210	TTGCTAAGCCACACCCCCAAGGGAT	ingroup-unique	-0.2804
   ...
   </pre>

2) FASTA file (.fa)

   A multi-FASTA with all sequences passing the filters set at run time.
   Input parameters will be captured in the filename.

   e.g.
   <pre>
    unikseq_v0-2-8beta-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-s100-p25-l1-u90-m0.fa

    -k length (option, default: -k 25)
    -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
    -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
    -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
    -u min. [% unique] kmers in regions (option, default: -u 90 %)
    -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m 0 % [original behaviour])
   </pre>

   In this example, unique sequences -s >=100 bp, found in -p >=25% of ingroup sequence entries on average, with -u >=90% of its -k 25-mers uniquely found (i.e. not in outgroup entries kmers), and with a leniency of at most -l 1 consecutive non-unique kmer (i.e. found in outgroup). Further, reference kmers are not tolerated in the outgroup (i.e. -m 0% outgroup entries are tolerated to have the reference kmer at each position).

   The header of each FASTA entry captures several key information.
   e.g.
   <pre>
   >KF597303.1region0-208_size209_propspcIN40.7_propunivsOUT99.5_avgOUTentries0.0
   GCTAGTGTAGCTTAATGTAAAG....CACGCACGTAGCCCAAGACAC

   1. The reference accession and start-end positions of the unique region (0-based coordinates, regionXX-XX)
   2. The size of the region (sizeXX) in bp
   3. The average proportion (%) of ingroup entries with reference kmers (over unique region length, propspcINXX.X %)
   4. The percent kmer uniqueness over the region length (propunivsOUTXX.X %)
   5. Average number of outgroup entries over the region length (avgOUTentriesXX.X). The lower the number the more unique the sequence.  
   </pre>   


### Plotting "butterfly" plots
-----------

Refer to example.r and replace these lines:

dfa<-read.table("XX unikseq-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-uniqueKmers.tsv XX", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("XX A. fragilis XX"), " Mt genome"))



