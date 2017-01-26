These are the perl scripts for making contigs from the 2 X 100 reads.
The directories are hard wired, so you will need to substitutue your own paths.

The input required is:
Unfiltered sequence reads (We used those generated from HLA homozygous individuals)
A file listing all possible 101mers from region of interest (in our case, the MHC project haplotypes )
A single reference sequence of the region for integrity checks 

Setup.pl calls the scripts that perform the following:

	selects reads for running in Velvet
 	runs Velvet
	Polishes ends
	Splits reads at homopolymer stretches
	Extends contigs
	runs Minimus Assembly
	
The output is a fasta file of contigs.