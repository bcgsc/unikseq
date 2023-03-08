[![Release](https://img.shields.io/github/release/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/unikseq/total?logo=github)](https://github.com/bcgsc/unikseq/releases/download/v1.0.0/unikseq-1.0.0.tar.gz)
[![Issues](https://img.shields.io/github/issues/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/issues)
Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/stargazers)

![Logo](https://github.com/bcgsc/unikseq/blob/main/unikseq-logo.png)

# unikseq
## Unique (& conserved) region identification in DNA sequences, using a k-mer approach
### 2020-2023


## Contents

1. [Description](#description)
2. [Implementation and requirements](#implementation)
3. [Install](#install)
4. [Dependencies](#dep)
5. [Bloom filter version](#bloom)
6. [Documentation](#docs)
7. [Citing unikseq](#cite)
8. [Credits](#credits)
9. [Running unikseq](#run)
10. [Test data](#data)
11. [Output](#output)
12. [Algorithm](#algorithm)
13. [Quick reference](#quickref)
14. [Generating "butterfly" plots](#bplot)
15. [License](#license)


## Description <a name=description></a>

Unikseq systematically processes the k-mers of a reference sequence, tolerated in an ingroup, but not (or marginally) tolerated in an outgroup sequence set to ultimately help identify regions that are unique within that reference. 

The unique sequences identified by unikseq are useful for designing qPCR primer-probe sets with high specificity, for instance, with no manual intervention nor time-consuming sequence inspection.

Unikseq has broad applications in genomics and biology, including species monitoring and conservation; It was used to develop successful environmental DNA (eDNA) mitogenome assays that are highly specific to their intended target, fast and efficiently.

Because unikseq does not rely on sequence alignments, it is much faster than multiple sequence alignments (MSA), and doesn't require additional downstream analyses one would need to carry out atfer MSA to identify unique, and potentially conserved, regions. Further, because unikseq employs a k-mer approach, in/outgroup FASTA sequence files need not be structured. As such, the input in/outgroup is very flexible and can include raw RNA-seq/WGS sequencing reads, contiguous/fragmented genome sequences with inconsistent start coordinates, unordered/unoriented contigs or even a mix of reads/genome sequences - as long as the genome sequences to compare a reference against are represented in full (i.e. complete), to the best of the user's knowledge; This is especially important for outgroup sequence sets, as absence of k-mers due to incomplete sequences may result in the identification of unique sequences in the reference sequence under scrutiny (see NOTES below). 

The unikseq-Bloom utility is a code variant for processing Gbp-scale genomes/sequencing data sets. Please note that the initial implementation requires pre-built Bloom filters data structures (generated with the writeBloom.pl utility, provided with unikseq). These are regular, and not counting Bloom filters; As such, k-mers are not counted, and their presence/absence alone are used to infer uniqueness.

```diff
! NOTE1: The unique (and/or hypervariable) reference regions identified by unikseq are RELATIVE to the sequences provided in the outgroup file. Please ensure COMPLETENESS* of outgroup sequences to help interpretability of results or for use in downstream applications
```

```diff
! *NOTE2: An incomplete outgroup sequence (or sequences) may be used; In that case, the reference regions identified by unikseq could help characterize missing stretches in those outgroup sequences, for instance
```


## Implementation and requirements <a name=implementation></a>

Unikseq is developed in PERL and runs on any system where PERL is installed.


## Install <a name=install></a>

Clone and enter the unikseq directory.
<pre>
git clone https://github.com/bcgsc/unikseq
cd unikseq
</pre>


## Dependencies <a name=dep></a>

If PERL is installed on your system, and you're interested in the original unikseq code*, you're good to go (no additional libraries needed, nor dependencies).

<pre>
* As opposed to the Bloom filter version, see below
</pre>

## Bloom filter version <a name=bloom></a>

If you are interested in running the Bloom filter version on your system and the pre-built libraries do not work on your system, please follow these instructions on how to recompile*:

https://github.com/bcgsc/LINKS/tree/d215339720f39c04537d859d7ea4e962c81b8d53#instructions-for-building-the-bloomfilter-perl-module

*note: You will need to download a older version of the LINKS long read genome scaffolder to do so (links v1.8.7) available here:
https://github.com/bcgsc/LINKS/releases/download/v1.8.7/links_v1-8-7.tar.gz


## Documentation <a name=docs></a>

Refer to the README.md file on how to install and run unikseq
Manuscript in preparation


## Citing unikseq <a name=cite></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/stargazers) and for using, developing and promoting this free software!

If you use unikseq in your research, please cite: TBD


## Credits <a name=credits></a>

unikseq (concept, algorithm design and prototype): Rene Warren


## Running unikseq <a name=run></a>

<pre>
Usage: ./unikseq.pl VERSION
-----input files-----
 -r reference FASTA (required)
 -i ingroup FASTA (required)
 -o outgroup FASTA (required)
-----k-mer uniqueness filters-----
 -k length (option, default: -k 25)
 -l [leniency] min. non-unique consecutive k-mers allowed in outgroup (option, default: -l 0)
 -m max. [% entries] in outgroup tolerated to have a reference k-mer (option, default: -m 0 % [original behaviour])
-----output filters-----
 -t print only first t bases in tsv output (option, default: -t [k])
 -c output conserved FASTA regions between reference and ingroup entries (option, -c 1==yes -c 0==no, [default, original unikseq behaviour])
 -s min. reference FASTA region [size] (bp) to output (option, default: -s 100 bp)
 -p min. [-c 0:region average /-c 1: per position] proportion of ingroup entries (option, default: -p 0 %)
 -u min. [% unique] k-mers in regions (option, default: -u 90 %)
 -v output tsv files (option, -v 1==yes -v 0==no [default])

or

Usage: ./unikseq-Bloom.pl VERSION
 -i ingroup Bloom filter (required)
 -o outgroup Bloom filter (required)
The rest of the options are the same.

Bloom filter data structures are built with the writeBloom.pl utility included in the "tools"  directory.
Both writeBloom.pl and unikseq-Bloom.pl have a library dependency (mentioned in the "Install" section above), with the pre-compiled library "BloomFilter.pm  BloomFilter.so" located in the "libexec" folder, "tools" directory. 


</pre>

Notes:
<pre>

 -r reference FASTA (required)
  reference FASTA (unique FASTA) analysis is done relative to it. Could be a multi-FASTA, but avoid sequences that are too long, especially when running in -c 1 mode (ideally < 10Mbp).

 -i ingroup FASTA (required)
  tolerated sequences. Used to find regions unique to a % (see -p option).
  In v1.2.1 onward, multi-FASTA entries are grouped by the first non-space identifier in the FASTA header. This is useful when querying k-mers from genome assemblies or even sequencing reads, summarizing for a giving species, for instance. e.g. >myID contig1 and >myID contig2. When unikseq calculates proportions reported in output files, it will summarize by "myID", counting multi-FASTA entries as one. The original behaviour is obtained by using distinct, non-space headers for each entries (e.g. >myID_contig1 and >myID_contig2 count as two entries).

 -o outgroup FASTA (required)
  outgroup to query k-mers against. Note that input reference and ingroup sequences will be automatically excluded from this set.
  In v1.2.1 onward, multi-FASTA entries are grouped by the first non-space identifier in the FASTA header. This is useful when querying k-mers from genome assemblies or even sequencing reads, summarizing for a giving species, for instance. e.g. >myID contig1 and >myID contig2. When unikseq calculates proportions reported in output files, it will summarize by "myID", counting multi-FASTA entries as one. The original behaviour is obtained by using distinct, non-space headers for each entries (e.g. >myID_contig1 and >myID_contig2 count as two entries). For those interested in sequence conservation between reference and ingroup only, you have the option to supply a dummy -o outgroup FASTA file consisting of a header and single base, for instance. 

 -k length (option, default: -k 25)
  k-mer length.

 -l [leniency] min. non-unique consecutive k-mers allowed in outgroup (option, default: -l 1)
  this leniency factor controls the tolerance for non-unique k-mers in the outgroup when computing a sequence stretch. If set to 0, the sequence stretch is limited to unique k-mers only and the % unique bases of the stretch will be 100% (see next parameter). Set to 1, allows some leniency and tolerates a single non-unique k-mer in a row. Set to 2, tolerates up to 2 in a row, etc. When set >0, the %unique bases of a sequence stretch will almost never be 100% unique and so -u MUST be lowered below 100%.

 -m max. [% entries] in outgroup tolerated to have reference k-mer (option, default: -m 0 % [original behaviour])
  controls the k-mer "uniqueness" in the outgroup, by tolerating a certain fraction of sequence entries in the outgroup having the reference k-mer. This option could be useful when there's high similarity between the reference, ingroup AND outgroup sequences and more fine-grain adjustments are needed. Not specifying this option (-m 0) is the original unikseq behaviour.

 -t print only first t bases in tsv output (option, default: -t [k])
  Avoids writing too much data to file. Users opt to specify the first -t base(s) to be printed in the tsv files, when the -v option is turned on (i.e. -v 1, see -v option below). Original (default) behaviour is to print to file the whole [k]-mer sequence.

 -c output conserved FASTA regions between reference and ingroup entries (option, -c 1==yes -c 0==no, [default, original unikseq behaviour])
  Boolean (1/0, yes/no), controls the output behaviour of unikseq. When set (-c 1), conserved regions between ingroup sequence entries and the reference are identified by repurposing the -p parameter (below), evaluating the % of entries having a reference k-mer at each position. When the % falls below the set -p threshold, regions will be part of a new unikseq FASTA output (.conserved.fa) and will be marked in upper-case (A,C,G,T) in the unique regions identified by unikseq. This option is useful for quick identification of unique sequences that are also conserved (to a tunable degree, and controlled by -p). The default (-c 0) is the original unikseq behaviour. Note: -c 1 may result in long run times on large (>10Mbp) reference (-r) sequences.

 -s min. reference FASTA region [size] (bp) to output (option, default: -s 100 bp)
  minimum "unique" (and -c 1:"conserved") reference (target) region size to report.

 -p min. [-c 0:region average /-c 1: per position] proportion of ingroup entries (option, default: -p 0 %)
  unikseq tracks the number of qualifying sequences in the ingroup over the sequence stretch, averages and calculates a proportion of the total entries in the ingroup. In the -c 0 mode (original unikseq behaviour), sequences are reported only when that proportion is above the minimum set. When -c 1 is set, -p does not impose a minimum average proportion of ingroup entries within unique regions and instead, non-conserved regions are soft-masked (a,c,g,t) and conserved regions are upper-cased (A,C,G,T) in the FASTA output. 

 -u min. [% unique] k-mers in regions (option, default: -u 90 %)
  controls for sequence uniqueness in the reference output regions. 

 -v output tsv files (option, -v 1==yes -v 0==no [default])
  Boolean (1/0, yes/no), controls tsv file output. By default, tsv files are no longer output by unikseq. To turn it on, set it to 1 (e.g. -v 1)


 Example command (refer to the Test data section below):
 ./unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o actinopterygii.fa -s 100 -p 25 -l 1 -u 90

</pre>

### Test data <a name=data></a>
---------
We showcase the utility of unikseq in identifying mitogenome (MtG) regions unique to C. maximum (basking shark, CEMA.fa reference) compared to ray-finned fishes (actinopterygii.fa outgroup, non-target sequence set n=868 MtG), and conserved in >=25% of (shark.fa ingroup, n=189 MtG) shark mitogenomes, on average.

Depending on your system, expect unikseq to use 2.5 GB RAM and run in 49.1s (wall clock time) on a single CPU thread on a MacBook Pro (2.6 GHz 6-Core Intel Core i7 chipset with 16GB RAM onboard) running mac OS (Catalina v10.15.7). On a server-class CentOS Linux 7 system with 144 Intel(R) Xeon(R) Gold 6254, 3.10GHz CPUs with 3TB RAM, the same test sample ran in 31.5s (wall clock time), using a single thread and required 2.5GB RAM. If your system is limited in RAM, you could subsample from shark.fa and/or actinopterygii.fa -- just to make sure unikseq is installed properly and will run on your system.

<pre>
1. Go to ./testdata
(cd testdata)

2. Unzip all FASTA files
(gunzip *fa) on unix

3. Run unikseq on the provided test data
../unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o actinopterygii.fa -s 100 -p 25 -l 1 -u 90

This specific command will generate two output files:
unikseq_vXX-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-c0-s100-p25-l1-u90-m0-unique.fa
unikseq_vXX-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-c0-s100-p25-l1-u90-m0-unique.log

If the run is successful, the -unique.fa FASTA output should contain 5 sequences.
</pre>


## Output (-c0 : 3 files / -c1 : 5 files) <a name=output></a>

1) TSV file (-uniqueKmers.tsv) 

   Tab-Separated Variable file. Reports all reference sequence k-mers in 4 columns:
   <pre>
   header	position	[k]-mer	condition	value
   [FASTA header][coordinates][sequence][in/out group][proportion in each in/out group]
   By default, every instance of a reference k-mer is reported when found in the outgroup.
   When it is not found, the ingroup-unique will be reported (if found). Note: when -t is specified, only the first -t bases of the k-mer will be shown in the tsv file(s), but the data reported is for the whole k-mer.
   If a reference k-mer is found in outgroup sequences, the ingroup-unique WILL NOT report any
   values. This is the file to use to generate "butterfly" plots (see below and r script attached with this distribution [example.r] for details)

   e.g.
   header  position        25-mer  condition       value
   KC914387.1	0	GCTAGTGTAGCTTAATGTAAAGTAT	ingroup-unique	-0.1270
   ...
   KC914387.1	207	ACCTTGCTAAGCCACACCCCCAAGG	ingroup-unique	-0.9683
   KC914387.1	208	CCTTGCTAAGCCACACCCCCAAGGG	outgroup	0.0115
   KC914387.1	209	CTTGCTAAGCCACACCCCCAAGGGA	outgroup	0.0115
   KC914387.1	210	TTGCTAAGCCACACCCCCAAGGGAT	ingroup-unique	-0.2804
   ...
   </pre>

2) FASTA file (-unique.fa)

   A multi-FASTA with all unique (relative to outgroup) sequences passing the filters set at run time.
   Input parameters will be captured in the filename.

   e.g.
   <pre>
    unikseq_v1.1.0beta-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-c0-s100-p25-l1-u90-m0.fa

    -k length (option, default: -k 25)
    -c conserved mode (option, default: -c 0)
    -s min. reference region [size] (bp) to output (option, default: -s 100 bp)
    -p min. average [proportion] ingroup entries in regions (option, default: -p 25 %)
    -l [leniency] min. non-unique consecutive k-mers allowed in outgroup (option, default: -l 1)
    -u min. [% unique] k-mers in regions (option, default: -u 90 %)
    -m max. [% entries] in outgroup tolerated to have reference k-mer at each position (option, default: -m 0 % [original behaviour])
   </pre>

   In this example, unique sequences -s >=100 bp, found in -p >=25% of ingroup sequence entries on average, with -u >=90% of its -k 25-mers uniquely found (i.e. not in outgroup entries k-mers), and with a leniency of at most -l 1 consecutive non-unique k-mer (i.e. found in outgroup). Further, reference k-mers are not tolerated in the outgroup (i.e. -m 0% outgroup entries are tolerated to have the reference k-mer at each position).

   The header of each FASTA entry captures several key information.
   e.g.
   <pre>
   >KF597303.1region0-208_size209_propspcIN40.7_propunivsOUT99.5_avgOUTentries0.0
   GCTAGTGTAGCTTAATGTAAAG....CACGCACGTAGCCCAAGACAC

   1. The reference accession and start-end positions of the unique region (0-based coordinates, regionXX-XX)
   2. The size of the region (sizeXX) in bp. It will be equal or larger than -s.
   3. The average proportion (%) of ingroup entries with reference k-mers (over unique region length, propspcINXX.X %). It will be equal or higher than -p.
   4. The percent k-mer uniqueness over the region length (propunivsOUTXX.X %). It will be equal or higher than -u.
   5. Average number of outgroup entries over the region length (avgOUTentriesXX.X). The lower the number the more unique the sequence.
   </pre>   

   ```diff
   ! NOTE: When -c 1 is set, -p does not impose a minimum average proportion of ingroup entries within unique regions and instead, non-conserved regions are soft-masked (a,c,g,t) and conserved regions are upper-cased (A,C,G,T) in the FASTA output. This handy feature enables quick identification of conserved (ingroup) regions within unique (relative to outgroup) sequences.
   ```

3) LOG file (.log)

   Captures the verbose STDOUT in a log file, showing the progress of the program.    

   ```diff
   ! NOTE: When -c is set (-c 1), there are two additional output files:
   ```

4) TSV file (-conservedKmers.tsv)

   Tab-Separated Variable file. Reports all reference sequence k-mers in 5 columns:
   <pre>
   header	position     [k]-mer    condition       num_entries	proportion
   [coordinates][sequence][ingroup][number of entries in ingroup with conserved reference k-mer][proportion relative to ingroup entries]. Note: when -t is specified, only the first -t bases of the k-mer will be shown in the tsv file(s), but the data reported is for the whole k-mer.

   e.g.
   header	position        25-mer    condition       num_entries     proportion
   ...
   KC914387.1	16      TTTAAAGCATGGCACTGAAGATGCT       ingroup 121     1.0000
   KC914387.1	17      TTAAAGCATGGCACTGAAGATGCTA       ingroup 121     1.0000
   KC914387.1	18      TAAAGCATGGCACTGAAGATGCTAA       ingroup 121     1.0000
   KC914387.1	19      AAAGCATGGCACTGAAGATGCTAAG       ingroup 115     0.9504
   KC914387.1	20      AAGCATGGCACTGAAGATGCTAAGA       ingroup 115     0.9504
   KC914387.1	21      AGCATGGCACTGAAGATGCTAAGAT       ingroup 115     0.9504
   ...
   </pre>


5) FASTA file (-conserved.fa)

   A multi-FASTA with all conserved (within ingroup) sequences passing the filters set at run time.
   Input parameters are captured in the filename (see above description in (2)).

   The header of each FASTA entry captures several key information.
   e.g.
   <pre>
   >ch-CACA1-CM030070.1region0-1527_size1528_propspcIN99.8
   ...

   1. The reference accession and start-end positions of the conserved region (0-based coordinates, regionXX-XX)
   2. The size of the region (sizeXX) in bp. It will be equal or larger than -s.
   3. The average proportion (%) of ingroup entries with reference k-mers (over unique region length, propspcINXX.X %). It will be equal or higher than -p.
   </pre>


## Algorithm design and implementation <a name=algorithm></a>

![unikseqMethod](https://github.com/bcgsc/unikseq/blob/main/unikseq-method.png)
The algorithm starts by first parsing FASTA sequence(s) supplied by the user as “outgroup” (-o option) and “ingroup” (-i option) and extracting every word of length k (k-mers, -k option) and their reverse complement and storing each in respective two-dimensional hash (or Bloom filter) data structures, keeping track of the k-mer occurrence in each FASTA entry for either sets. We point out that in/outgroup sequences need not start at the same position, nor be represented on the same strand since unikseq is k-mer based, so no specific DNA sequence formatting is required other than supplying a FASTA-formatted file (e.g. no need to adjust the sequence start for mtDNA genomes). 

The search for stretches of potentially unique sequence in a reference sequence (-r option) begins with the 5’ to 3’ extraction of a forward-strand k-mer and querying of the above two data structures for 1) exclusion of the k-mer in the outgroup and 2) inclusion of the k-mer in the ingroup, moving the k-mer frame over base by base until the entire FASTA sequence is read and all reference k-mers have been interrogated. When a sequence k-mer is not found in the outgroup data structure, it is deemed unique and its position and coverage in the outgroup set is tracked and a unique sequence stretch is initialized if the k-mer demarks the beginning of a new unique region. A record of the k-mer’s presence in the ingroup data structure is also kept for the purpose of plotting and evaluating the overall conservation of the unique sequence stretch in the ingroup (see the -p option described below). In the event that a unique k-mer follows a k-mer previously identified as unique, the last unique base is concatenated onto the growing unique sequence stretch on the 3’-end. This process is repeated until the end of the FASTA reference sequence is reached or until a condition is no longer met, including when the next k-mer is found to be non-unique in the outgroup set. 

To facilitate the detection of unique regions that may be interspersed with non-unique k-mers, a leniency parameter (-l option) is used to control the number of consecutive non-unique reference k-mers that are tolerated in the outgroup. The unique sequence stretch extension is stopped when the number of consecutive reference k-mers found in the outgroup set has exceeded that minimum -l threshold. A grace parameter (-m option) can also be used to adjust the tolerance, or uniqueness factor, in the outgroup. This is particularly useful when searching for unique regions in sets of very similar sequences, for instance, as -m describes the maximum proportion of outgroup sequences in the set that are tolerated before the reference k-mer is considered unique (default set to 0, indicating no tolerance). The -m and -l parameters work hand-in-hand in controlling the stringency of unikseq, with -m describing the unique k-mer tolerance rate threshold in the outgroup set, and the latter option allowing -l consecutive “non-unique” bases when that is not the case.

When the unique sequence can no longer be extended, it will be written to a FASTA file only when 1) it has reached a length of at least (-s option) bp in size, 2) its constituent k-mers have been identified, on average proportion, in at least (-p option) % of ingroup sequences and, overall, 3) those k-mers are at least (-u option) % unique, as specified by -l and -m. 

Unikseq also outputs a tab-separated value (tsv) file that tracks, at each coordinate relative to the reference sequence, the proportion of corresponding k-mers in the outgroup and ingroup sets, and was used to generate the butterfly plots (below). Instructions and code for generating the butterfly plots in the R programming language are both available from the repository at the URL below. The intended use-case of unikseq is for identification of unique sequences in mitochrondrial [genome] and similarly short sequences (In initial tests comparing 52 x 2-5Mbp bacterial [HMP mock community] genomes, unikseq ran in 5m and required 54GB RAM). For larger, Gbp-size genomes, we recommend the use of unikseq-Bloom, included in the v1.3 commit of the code. Unikseq is developed in PERL and runs on any system where PERL is installed, requiring no additional library. It is distributed under GPLv3 license and available freely from github (https://github.com/bcgsc/unikseq).

 
## Quick reference <a name=quickref></a>

Below is a reference guide for controlling the [stringency &] output of unikseq.
<pre>
(+) more stringent (identifying less & potentially smaller regions, with increased uniqueness)

   (-) decrease: -l and -m (controls region uniqueness) parameters
   (+) increase: -p, -u and -s (output filter controls) parameters. You may also increase k

(-) less stringent (identifying more & potentially larger regions, with decreased uniqueness)

   (+) increase: -l and -m (controls region uniqueness) parameters
   (-) decrease: -p, -u and -s (output filter controls) parameters. You may also decrease k
</pre>

## Generating "butterfly" plots <a name=bplot></a>

![UnikseqButterflyPlot](https://github.com/bcgsc/unikseq/blob/main/unikseq-butterfly.png)

Refer to `example.r` included with the unikseq distribution, and replace these lines:

<pre>
dfa<-read.table("XX unikseq-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-uniqueKmers.tsv XX", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("XX C. maximus XX"), " Mt genome"))
</pre>

## License <a name=license></a>

Unikseq Copyright (c) 2020-2023 British Columbia Cancer Agency Branch.  All rights reserved.

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
