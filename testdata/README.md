23 February 2023

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
../unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90
---or---
../unikseq.pl -k 25 -r CEMA.fa.gz -i shark.fa.gz -o teleost.fa.gz -s 100 -p 25 -l 1 -u 90
</pre>
