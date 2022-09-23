[![Release](https://img.shields.io/github/release/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/releases)
[![Issues](https://img.shields.io/github/issues/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/issues)
Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/stargazers)

![Logo](https://github.com/bcgsc/unikseq/blob/main/unikseq-logo.png)

# unikseq
Unique (& conserved) region identification in DNA sequences, using a kmer approach

## Rene L. Warren, 2020-2022
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

Unikseq systematically processes the kmers of a reference sequence, tolerated in an ingroup, but not (or marginally tolerated) in an outgroup to ultimately help identify regions that are unique within that reference. 

The unique sequences identified by unikseq are useful for designing qPCR primer-probe sets with high specificity, for instance, with no manual intervention nor sequence inspection.

Unikseq has broad applications in genomics and biology, and was used to develop successful environmental DNA (eDNA) mitogenome assays that are highly specific to their intended target, fast.

Because unikseq does not rely on sequence alignments, it is much faster than multiple sequence alignments (MSA), and doesn't require the additional analysis one would need to carry out atfer an MSA to identify unique conserved regions.

Unikseq is developed in PERL and runs on any system where PERL is installed.


### Running unikseq
-----------

<pre>
Usage: ./unikseq.pl v1.0.0
-----input files-----
 -r reference FASTA (required)
 -i ingroup FASTA (required)
 -o outgroup FASTA (required)
-----kmer uniqueness filters-----
 -k length (option, default: -k 25)
 -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
 -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m 0 % [original behaviour])
-----output filters-----
 -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
 -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
 -u min. [% unique] kmers in regions (option, default: -u 90 %)
</pre>

Notes:
<pre>

 -r reference FASTA (required)
  reference FASTA (unique FASTA) analysis is done relative to it

 -i ingroup FASTA (required)
  tolerated sequences. Used to find regions unique to a % (see -p option)

 -o outgroup FASTA (required)
  outgroup to query kmers against. Note that input reference and ingroup sequences will be automatically excluded from this set.

 -k length (option, default: -k 25)
  kmer length

 -l [leniency] min. non-unique consecutive kmers allowed in outgroup (option, default: -l 1)
  this leniency factor controls the tolerance for non-unique kmers in the outgroup when computing a sequence stretch. If set to 0, the sequence stretch is limited to unique kmers only and the % unique bases of the stretch will be 100% (see next parameter). Set to 1, allows some leniency and tolerates a single non-unique kmer in a row. Set to 2, tolerates up to 2 in a row, etc. When set >0, the %unique bases of a sequence stretch will almost never be 100% unique and so -u MUST be lowered below 100%.

 -m max. [% entries] in outgroup tolerated to have reference kmer at each position (option, default: -m 0 % [original behaviour])
  controls the kmer "uniqueness" in the outgroup, by tolerating a certain fraction of sequence entries in the outgroup having the reference kmer. This option could be useful when there's high similarity between the reference, ingroup AND outgroup sequences and more fine-grain adjustments are needed. Not specifying this option (-m 0) is the original unikseq behaviour.

 -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
  minimum "unique" reference (target) region size to report.

 -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
  program tracks the number of qualifying sequences in the ingroup over the sequence stretch, averages and calculates a proportion of the total entries in the ingroup. Sequences are reported only when that proportion is above the minimum set here. 

 -u min. [% unique] kmers in regions (option, default: -u 90 %)
  controls for sequence uniqueness in the reference output regions. 

 Example command:
 ./unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90

</pre>


### Outputs (3)

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
   position        kmer                 condition       value
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
    unikseq_v1.0.0-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-s100-p25-l1-u90-m0.fa

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
   2. The size of the region (sizeXX) in bp. It will be equal or larger than -s.
   3. The average proportion (%) of ingroup entries with reference kmers (over unique region length, propspcINXX.X %). It will be equal or higher than -p.
   4. The percent kmer uniqueness over the region length (propunivsOUTXX.X %). It will be equal or higher than -u.
   5. Average number of outgroup entries over the region length (avgOUTentriesXX.X). The lower the number the more unique the sequence.
   </pre>   

3. LOG file (.log)

   Captures the verbose STDOUT in a log file, showing the progress of the program.    


### Unikseq algorithm design and implementation

The algorithm starts by first parsing FASTA sequence(s) supplied by the user as “outgroup” (-o option) and “ingroup” (-i option) and extracting every word of length k (kmers, -k option) and their reverse complement and storing each in respective two-dimensional hash data structures, keeping track of the kmer occurrence in each FASTA entry for either sets. We point out that in/outgroup sequences need not start at the same position, nor be represented on the same strand since unikseq is kmer-based, so no specific DNA sequence formatting is required other than supplying a FASTA-formatted file (e.g. no need to adjust the sequence start for mtDNA genomes). 

The search for stretches of potentially unique sequence in a reference sequence (-r option) begins with the 5’ to 3’ extraction of a forward-strand kmer and querying of the above two hash data structures for 1) exclusion of the kmer in the outgroup and 2) inclusion of the kmer in the ingroup, moving the kmer frame over base by base until the entire FASTA sequence is read and all reference kmers have been interrogated. When a sequence kmer is not found in the outgroup hash data structure, it is deemed unique and its position and coverage in the outgroup set is tracked and a unique sequence stretch is initialized if the kmer demarks the beginning of a new unique region. A record of the kmer’s presence in the ingroup hash data structure is also kept for the purpose of plotting and evaluating the overall conservation of the unique sequence stretch in the ingroup (see the -p option described below). In the event that a unique kmer follows a kmer previously identified as unique, the last unique base is concatenated onto the growing unique sequence stretch on the 3’-end. This process is repeated until the end of the FASTA reference sequence is reached or until a condition is no longer met, including when the next kmer is found to be non-unique in the outgroup set. 

To facilitate the detection of unique regions that may be interspersed with non-unique kmers, a leniency parameter (-l option) is used to control the number of consecutive non-unique reference kmers that are tolerated in the outgroup. The unique sequence stretch extension is stopped when the number of consecutive reference kmers found in the outgroup set has exceeded that minimum -l threshold. A grace parameter (-m option) can also be used to adjust the tolerance, or uniqueness factor, in the outgroup. This is particularly useful when searching for unique regions in sets of very similar sequences, for instance, as -m describes the maximum proportion of outgroup sequences in the set that are tolerated before the reference kmer is considered unique (default set to 0, indicating no tolerance). The -m and -l parameters work hand-in-hand in controlling the stringency of unikseq, with -m describing the unique kmer tolerance rate threshold in the outgroup set, and the latter option allowing -l consecutive “non-unique” bases when that is not the case.

When the unique sequence can no longer be extended, it will be written to a FASTA file only when 1) it has reached a length of at least (-s option) bp in size, 2) its constituent kmers have been identified, on average proportion, in at least (-p option) % of ingroup sequences and, overall, 3) those kmers are at least (-u option) % unique, as specified by -l and -m. 

Unikseq also outputs a tab-separated value (tsv) file that tracks, at each coordinate relative to the reference sequence, the proportion of corresponding kmers in the outgroup and ingroup sets, and was used to generate the butterfly plots (below). Instructions and code for generating the butterfly plots in the R programming language are both available from the repository at the URL below. The intended use-case of unikseq is for identification of unique sequences in mitochrondrial [genome] or shorter sequences, and thus has not been tested on larger genomes. Unikseq is developed in PERL and runs on any system where PERL is installed, requiring no additional library. It is distributed under GPLv3 license and available freely from github (https://github.com/bcgsc/unikseq).

 
### Quick reference

Below is a reference guide for controlling the [stringency &] output of unikseq.
<pre>
(+) more stringent (identifying less & potentially smaller regions, with increased uniqueness)

   (-) decrease: -l and -m (controls region uniqueness) parameters
   (+) increase: -p, -u and -s (output filter controls) parameters. You may also increase k

(-) less stringent (identifying more & potentially larger regions, with decreased uniqueness)

   (+) increase: -l and -m (controls region uniqueness) parameters
   (-) decrease: -p, -u and -s (output filter controls) parameters. You may also decrease k
</pre>

### Plotting "butterfly" plots
-----------

![UnikseqButterflyPlot](https://github.com/bcgsc/unikseq/blob/main/unikseq-butterfly.png)

Refer to `example.r` included with the unikseq distribution, and replace these lines:

<pre>
dfa<-read.table("XX unikseq-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-uniqueKmers.tsv XX", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("XX A. fragilis XX"), " Mt genome"))
</pre>

### License
-------

Unikseq Copyright (c) 2020-2022 British Columbia Cancer Agency Branch.  All rights reserved.

Unikseq is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
