# projectH

These are the scripts and accessory files from the MHC region sequencing project as described in:

Sequences of 95 human MHC haplotypes reveal extreme coding variation in genes other than highly polymorphic HLA class I and II.
Paul J. Norman, Steven J. Norberg, Lisbeth A. Guethlein, Neda Nemat-Gorgani, Thomas Royce, Emily E. Wroblewski, Tamsen Dunn, Tobias Mann, Claudia Alicata, Jill A. Hollenbach, Weihua Chang, Melissa Shults Won, Kevin L. Gunderson, Laurent Abi-Rached, Mostafa Ronaghi and Peter Parham.
Genome Research 2017

All sequencing and assembled MHC haplotype data are linked to the NCBI BioProject (https://www.ncbi.nlm.nih.gov/bioproject/) accession number PRJEB6763

Included here are:

1. the scripts [see the scripts branch] and filters [filters branch] used to extract the CDS alleles for the expressed MHC region genes
	a. directly from fastq files (ExtractXX.bat)
	b. from the assembled contigs (Contigs_XX.bat)
and [data branch] the CDS alleles (also deposited in GenBank)

2. the scripts used to make contigs for each homozygous cell (Contig_generation/x.pl)
and [unprocessed_contigs] the resulting contigs

3. (Contigula) the scripts used to make full MHC region haplotype scaffolds from these contigs

4. Full alignments for every gene in the MHC region [alignments branch]
