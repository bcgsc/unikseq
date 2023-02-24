23 February 2023

Directory testdata

For convenience, we provide users with test data for running unikseq. The test data consists of:

reference (-r)
CEMA.fa
1 FASTA entry

ingroup sequence set (-i)
shark.fa
189 FASTA entries

outgroup/non-target sequence set (-o)
teleost.fa
872 FASTA entries

<pre>
1. Go to ./testdata
(cd testdata)

2. Unzip all FASTA files
(gunzip *fa) on unix

3. Run unikseq on the provided test data
../unikseq.pl -k 25 -r CEMA.fa -i shark.fa -o teleost.fa -s 100 -p 25 -l 1 -u 90
</pre>
