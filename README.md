[![Release](https://img.shields.io/github/release/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/unikseq/total?logo=github)](https://github.com/bcgsc/unikseq/releases/download/v2.0.1/unikseq-2.0.1.tar.gz)
[![Conda](https://img.shields.io/conda/dn/bioconda/unikseq?label=Conda)](https://anaconda.org/bioconda/unikseq)
[![Issues](https://img.shields.io/github/issues/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/issues)
[![link](https://img.shields.io/badge/unikseq-manuscript-brightgreen)](https://doi.org/10.1002/edn3.438)
[![link](https://img.shields.io/badge/unikseq2-manuscript-brightgreen)](https://doi.org/10.1111/1755-0998.70014)
Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/stargazers)

![Logo](https://github.com/bcgsc/unikseq/blob/main/unikseq-logo.png)

# unikseq
## Unique (& conserved) region identification in DNA sequences, using k-mers
### 2020-present


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

Unikseq has broad applications in genomics and biology, including species monitoring and conservation; It was [used to develop successful environmental DNA (eDNA) mitogenome assays that are highly specific to their intended target](https://doi.org/10.1002/edn3.438), fast and efficiently.

Because unikseq does not rely on sequence alignments, it is much faster than multiple sequence alignments (MSA), and doesn't require additional downstream analyses one would need to carry out atfer MSA to identify unique, and potentially conserved, regions. Further, because unikseq employs a k-mer approach, in/outgroup FASTA sequence files need not be structured. As such, the input in/outgroup is flexible and can include raw RNA-seq/WGS sequencing reads, contiguous/fragmented genome sequences with inconsistent start coordinates, unordered/unoriented contigs or even a mix of reads/genome sequences - as long as the genome sequences to compare a reference against are represented in full (i.e. complete), to the best of the user's knowledge; This is especially important for outgroup sequence sets, as absence of k-mers due to incomplete sequences may result in the identification of unique sequences in the reference sequence under scrutiny (see NOTES below). 

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

or using conda:
<pre>
conda install -c bioconda unikseq
</pre>


## Dependencies <a name=dep></a>

If PERL is installed on your system, and you're interested in the original unikseq code*, you're good to go (no additional libraries needed, nor dependencies).

<pre>
* As opposed to the Bloom filter version, see below
</pre>

## Bloom filter version <a name=bloom></a>

If you are interested in running the Bloom filter version on your system and the pre-built libraries* do not work on your system, please follow these instructions on how to recompile*:

https://github.com/bcgsc/LINKS/tree/d215339720f39c04537d859d7ea4e962c81b8d53#instructions-for-building-the-bloomfilter-perl-module

*note: You will need to download a older version of the LINKS long read genome scaffolder to do so (links v1.8.7) available here:
https://github.com/bcgsc/LINKS/releases/download/v1.8.7/links_v1-8-7.tar.gz

```diff
! *NOTE: The aforementioned conda install will install all dependencies, including the required Bloom filter libraries. We therefore highly recommend installing unikseq with conda (conda install -c bioconda unikseq).
```


## Documentation <a name=docs></a>

Refer to the README.md file on how to install and run unikseq. Read the 
[![link](https://img.shields.io/badge/unikseq-manuscript-brightgreen)](https://doi.org/10.1002/edn3.438) for additional implementation details.
![poster](https://github.com/bcgsc/unikseq/blob/main/unikseq-recomb2023poster.png)
Warren RL, Allison MJ, Lopez ML, Acharya-Patel N, Coombe L, Yang CL, Helbing CC, Birol I, "Unique region identification in genomes using a k-mer approach", 27th Annual International Conference on Research in Computational Molecular Biology (RECOMB), RECOMB, RECOMB-Seq and RECOMB-CG, Istanbul, Turkey, April 2023.


## Citing unikseq <a name=cite></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/unikseq.svg)](https://github.com/bcgsc/unikseq/stargazers) and for using, developing and promoting this free software!

If you use unikseq in your research, please cite: 
(1)
[Enabling robust environmental DNA assay design with “unikseq” for the identification of taxon-specific regions within whole mitochondrial genomes](https://doi.org/10.1002/edn3.438)
<pre>
Allison MJ, Warren RL, Lopez ML, Acharya-Patel N, Imbery JJ, Coombe L, Yang CL, Birol I, Helbing CC.
Enabling robust environmental DNA assay design with “unikseq” for the identification of 
taxon-specific regions within whole mitochondrial genomes. 
Environmental DNA, 00, 1-16 (2023)
https://doi.org/10.1002/edn3.438
</pre>
[![link](https://img.shields.io/badge/unikseq-manuscript-brightgreen)](https://doi.org/10.1002/edn3.438)

AND
(2)
[Conserved sequence identification within large genomic datasets using unikseq2: Application in environmental DNA assay development](https://doi.org/10.1111/1755-0998.70014)
<pre>
Lopez MLD, Warren RL, Allison MJ, Coombe L, Imbery JJ, Birol I, Helbing CC.
Conserved sequence identification within large genomic datasets using unikseq2: Application in environmental DNA assay development. Molecular Ecology Resources. 
Molecular Ecology Resources (2025)
https://doi.org/10.1111/1755-0998.70014
</pre>
[![link](https://img.shields.io/badge/unikseq2-manuscript-brightgreen)](https://doi.org/10.1111/1755-0998.70014)

## Credits <a name=credits></a>
<pre>
Unikseq concept, algorithm design, implementation: Rene L Warren
Project inception (iTrackDNA): Caren C Helbing
Project oversight (bioinformatics): Inanc Birol
Beta testers, eDNA assay designs: Michael J Allison, M. Louie Lopez, Neha Acharya-Patel
</pre>


## Running unikseq <a name=run></a>

<pre>
Usage: ./unikseq.pl
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

Usage: ./unikseq-Bloom.pl
 -i ingroup Bloom filter (required)
 -o outgroup Bloom filter (required)
The rest of the options are the same.

Bloom filter data structures are built with the writeBloom.pl utility included in the "tools"  directory.
Both writeBloom.pl and unikseq-Bloom.pl have a library dependency (mentioned in the "Install" section above), with the pre-compiled library "BloomFilter.pm  BloomFilter.so" located in the "libexec" folder, "tools" directory. 


</pre>

Notes:
<pre>

 -r reference FASTA (required)
  reference FASTA (unique FASTA) analysis is done relative to it. Could be a multi-FASTA, but avoid sequences that are too long, especially when running in -c 1 mode (ideally < 10Mbp). For longer sequences, use unikseq-Bloom.pl located in the tools folder. Version 1.3.4 supports gzip and zip files.

 -i ingroup FASTA (required)
  tolerated sequences. Used to find regions unique to a % (see -p option). Version 1.3.4 supports gzip and zip files.
  In v1.2.1 onward, multi-FASTA entries are grouped by the first non-space identifier in the FASTA header. This is useful when querying k-mers from genome assemblies or even sequencing reads, summarizing for a giving species, for instance. e.g. >myID contig1 and >myID contig2. When unikseq calculates proportions reported in output files, it will summarize by "myID", counting multi-FASTA entries as one. The original behaviour is obtained by using distinct, non-space headers for each entries (e.g. >myID_contig1 and >myID_contig2 count as two entries).

 -o outgroup FASTA (required)
  outgroup to query k-mers against. Note that input reference and ingroup sequences will be automatically excluded from this set. Version 1.3.4 supports gzip and zip files.
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
 ./unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90
 ---or---
 ./unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 25 -l 1 -u 90
</pre>

### Test data <a name=data></a>
---------
We showcase the utility of unikseq in identifying mitogenome (MtG) regions unique to C. maximum (basking shark, CEMA.fa reference) compared to ray-finned fishes (teleost.fa outgroup, non-target sequence set n=872 MtG), and conserved in >=25% of (shark.fa ingroup, n=189 MtG) shark mitogenomes, on average (unikseq.pl query below).

Depending on your system, expect unikseq.pl to use 2.5 GB RAM and run in 49.1s (wall clock time) on a single CPU thread on a MacBook Pro (2.6 GHz 6-Core Intel Core i7 chipset with 16GB RAM onboard) running mac OS (Catalina v10.15.7). On a server-class CentOS Linux 7 system with 144 Intel(R) Xeon(R) Gold 6254, 3.10GHz CPUs with 3TB RAM, the same test sample ran in 31.5s (wall clock time), using a single thread and required 2.5GB RAM. If your system is limited in RAM, you could subsample from shark.fa and/or teleost.fa -- just to make sure unikseq is installed properly and will run on your system.

#### Testing unikseq.pl
<pre>

1. Go to ./testdata
(cd testdata)


2. Unzip all FASTA files
(gunzip *fa.gz) on unix

Note: In version 1.3.4 and subsequent, unikseq.pl supports .zip and .gz FASTA files directly (i.e., no need to uncompressed before running) 


3. Run unikseq on the provided test data
../unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90
---or---
./unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 25 -l 1 -u 90 -v 1


This specific command will generate 3 output files (4 if you ran with -v 1):
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-c0-s100-p25-l1-u90-m0-unique.fa
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-c0-s100-p25-l1-u90-m0.bed
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-c0-s100-p25-l1-u90-m0.log
---or---
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-c0-s100-p25-l1-u90-m0-unique.fa
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-c0-s100-p25-l1-u90-m0.bed
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-c0-s100-p25-l1-u90-m0.log
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-uniqueKmers.tsv

If the run is successful, the -unique.fa FASTA output should contain 5 sequences.


4. Generate a butterfly plot with R (you may need to install rscript and the R package ggplot2)
rscript ../butterfly-plot.r unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-uniqueKmers.tsv "C. maximus"

5. Run unikseq in conserved (-c 1) mode on the provided test data
../unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 15 -l 1 -u 90 -v 1 -c 1
*note: this will output conserved regions (-conserved.fa) that satisfy both the minimum length (-s 100bp and up) threshold, and the proportion of conserved k-mers (-p 15% and higher) out of 189 shark entries in shark.fa.gz. The unique regions (-unique.fa) lists all unique regions in reference CEMA.fa.gz (i.e., not found in the 872 teleost.fa.gz outgroup entries). The conserved regions with shark.fa.gz (at the 15% level and higher) are represented as upper case bases in the latter file.

</pre>

#### Testing unikseq-Bloom.pl
<pre>

1. Go to ./testdata
(cd testdata)


2. Unzip all FASTA files
(gunzip *.fa.gz) on unix

Note: In version 1.3.4 and subsequent, unikseq-Bloom.pl supports .zip and .gz FASTA files directly (i.e., no need to uncompressed before running)


3. Generate Bloom filters from in/outgroup FASTA files

writeBloom.pl -f shark.fa
writeBloom.pl -f teleost.fa
---or---
writeBloom.pl -f shark.fa.gz
writeBloom.pl -f teleost.fa.gz


These commands will write two output files:
shark.fa_k25_p0.001_rolling.bloom
teleost.fa_k25_p0.001_rolling.bloom


4. Run unikseq-Bloom.pl on the provided test data and Bloom filters generated in (3)

unikseq-Bloom.pl -k 25 -r CEMA.fa -i shark.fa_k25_p0.001_rolling.bloom -o teleost.fa_k25_p0.001_rolling.bloom -s 500 -p 100 -l 1 -u 99.5
---or---
unikseq-Bloom.pl -k 25 -r CEMA.fa.gz -i shark.fa_k25_p0.001_rolling.bloom -o teleost.fa_k25_p0.001_rolling.bloom -s 500 -p 100 -l 1 -u 99.5


This specific command will generate three output files:
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0-unique.fa
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0.bed
unikseq_v2.0.1-r_CEMA.fa-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0.log
---or---
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0-unique.fa
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0.bed
unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa_k25_p0.001_rolling.bloom-o_teleost.fa_k25_p0.001_rolling.bloom-k25-c0-s500-p100-l1-u99.5-m0.log

If the run is successful, the -unique.fa FASTA output should contain 4 sequences.

</pre>

```diff
! NOTE: The sequence output of unikseq.pl and unikseq-Bloom.pl are NOT EXPECTED to be the same. This is because k-mers are counted (and counts used for analysis) in the former, but not in the latter
```


## Output (-c0 -v1 : 4 files / -c1 -v1 : 6 files) <a name=output></a>

1) TSV file (-uniqueKmers.tsv) 

   Tab-Separated Variable file. Reports all reference sequence k-mers in 5 columns:
   <pre>
   header	position	[k]-mer	condition	value
   [FASTA header][coordinates][sequence][in/out group][proportion in each in/out group]
   By default, every instance of a reference k-mer is reported when found in the outgroup.
   When it is not found, the ingroup-unique will be reported (if found). Note: when -t is specified, only the first -t bases of the k-mer will be shown in the tsv file(s), but the data reported is for the whole k-mer.
   If a reference k-mer is found in outgroup sequences, the ingroup-unique WILL NOT report any
   values. This is the file to use to generate "butterfly" plots (see below and r script attached with this distribution [butterfly-plot.r] for details)

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
    unikseq_v1.1.0beta-r_CEMA.fa-i_shark.fa-o_teleost.fa-k25-c0-s100-p25-l1-u90-m0.fa

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

3) BED-like file (.bed)

   The first three columns contain the names of chromosomes or scaffolds, the start, and the end coordinates of the sequences considered. When -c 1 is set, the bed-like file records the coordinates of conserved (upper-case base) regions ONLY, when applicable. The file can be used in conjunction with bedtools (intersect) to find possible overlap with the annotated regions of a GFF3 file, for instance. Example command:

   ```diff
   bedtools intersect -a KF597303.1.gff3 -b unikseq_v1.3.5-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-c0-s100-p25-l1-u90-m0.bed
   ```

4) LOG file (.log)

   Captures the verbose STDOUT in a log file, showing the progress of the program.    

   ```diff
   ! NOTE: When -c is set (-c 1), there are two additional output files:
   ```

5) TSV file (-conservedKmers.tsv)

   Tab-Separated Variable file. Reports all reference sequence k-mers in 6 columns:
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

6) FASTA file (-conserved.fa)

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
For user convenience, an R script is included to facilitate data visualization. Refer to `butterfly-plot.r` included with the unikseq distribution. You may test it using the example data provided; Please refer to the testdata folder README.md for instructions on how to generate butterfly plots.
![UnikseqButterflyPlot](https://github.com/bcgsc/unikseq/blob/main/unikseq-butterfly.png)
Example butterfly-type plot generated from the testdata (and run parameters) provided.


## License <a name=license></a>

Unikseq Copyright (c) 2020-present British Columbia Cancer Agency Branch.  All rights reserved.

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
