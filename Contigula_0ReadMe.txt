Here are the three perl scripts used to build the assemblies from the contigs that were generated from HLA homozygous individuals.
The inputs are the contigs and the sequence reads (with  those mapping outside MHC removed)
The order of use is Reorder_pre, Condense, ScaffMaker

1. Sort the contigs to remove non-HLA region (e.g. use F001.bat)
2. Remove read pairs having both ends that map within a single contig (e.g. use F002.bat)
3. Harvest read pairs having each end in a different contig (e.g. use F003.bat)
4. Run Reorder_pre.pl
5. Run Condense.pl
6. Run ScaffMaker.pl
7. Run GapFiller 
8. Order scaffolds relative to reference (SuperCoxMix), and remove those having both ends map within another (unless the ends have been clipped to do this -as marked with S flag in CIGAR)
	(e.g. F004.bat and F005.bat)
9. Run GapFiller

GapFiller parameters used: -m 30 -o 2 -r 0.7 -n 18 -d 550 -t 10 -g 0 -T 5 -i 5
