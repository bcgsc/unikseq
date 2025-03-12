23 February 2023 (updated 12 March 2025)

Directory testdata

For convenience, we provide users with mitogenome test data for running unikseq. The test data consists of:

reference mitogenome (-r)
CEMA.fa
1 FASTA entry

ingroup mitogenome (tolerated) sequence set (-i)
shark.fa
189 FASTA entries

outgroup/non-target (non-tolerated) sequence set (-o)
teleost.fa
872 FASTA entries

<pre>
1. Go to ./testdata
(cd testdata)


2. Unzip all FASTA files
(gunzip *.fa.gz) on unix

Note: In version 1.3.4 and subsequent, unikseq.pl and unikseq-Bloom.pl support .zip and .gz FASTA files directly (i.e., no need to uncompressed before running)


3. Run unikseq on the provided test data
../unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90 -v 1
---or---
../unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 25 -l 1 -u 90 -v 1


4. Generate a butterfly plot with R (must have rscript and ggplot2 installed)
rscript ../butterfly-plot.r unikseq_v2.0.1-r_CEMA.fa.gz-i_shark.fa.gz-o_teleost.fa.gz-k25-uniqueKmers.tsv "unikseqTest"


5. Run unikseq in conserved (-c 1) mode on the provided test data
../unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 15 -l 1 -u 90 -v 1 -c 1

*note: this will output conserved regions (-conserved.fa) that satisfy both the minimum length (-s 100bp and up) threshold, and the proportion of conserved k-mers (-p 15% and higher) out of 189 shark entries in shark.fa.gz. The unique regions (-unique.fa) lists all unique regions in reference CEMA.fa.gz (i.e., not found in the 872 teleost.fa.gz outgroup entries). The conserved regions with shark.fa.gz (at the 15% level and higher) are represented as upper case bases in the latter file.

</pre>
