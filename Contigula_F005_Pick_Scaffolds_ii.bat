
	## makes final tiling of scaffs prior to GF2 step##
		
		##This needs e.g. list/COX.txt ##
		
	##also needs gapfilled.final.fa converted to linear .fas and placed in /scaffss ##
	
ls list/*.txt |sed 's|list/||' | sed 's|.txt||' >contigs.txt

REF=superCOXmix
echo NNNNNNNNNNNNNNNNNNNN > Ns.txt

while read CONTIGS
do
echo $CONTIGS
#CONTIGS=DBB #tester

rm -f test*

#cat -v scaffs/$CONTIGS* >test.in
tr -d '\015' < scaffs/$CONTIGS* > test.in

cat list/$CONTIGS".txt" >testscaffs.txt

	while read SCAFF
	do
echo $SCAFF

grep -w $SCAFF test.in > test.arse
paste test.arse Ns.txt > test.scaff

cat test.oot test.scaff > temp.arse
mv temp.arse test.oot

	done < testscaffs.txt

cut -f2,3 test.oot |echo $(cat)|sed 's| ||g' >test.seq.oot

echo ">"$CONTIGS"_on_"$REF >test.head
cat test.head test.seq.oot > doneOrdered/$CONTIGS"_on_"$REF".fas"

done < contigs.txt

rm -f ref*
rm -f *.arse
