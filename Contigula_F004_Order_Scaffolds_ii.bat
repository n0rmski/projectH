			## orders scaffs prior to GF2 step##
		
		##This is for e.g. contigs/COX_150_GF.gapfilled.final.fa ##

	##also needs gapfilled.final.fa converted to linear .fas and placed in /scaffss ##

ls contigs/*.fa |sed 's|contigs/||' | sed 's|.fa||' >contigs.txt

ls ref/*.fas |sed 's|ref/||' | sed 's|.fas||' >ref.txt

REF=superCOXmix

cp ref/$REF".fas" ref.fna
bwa index ref.fna

while read CONTIGS
do
echo $CONTIGS

rm -f test*

cat contigs/$CONTIGS".fa" > miniplus.fas

k=64
bwa mem -k$k ref.fna miniplus.fas -M >oot.sam

awk '$2 == 4' oot.sam > All_done/$CONTIGS"vs"$REF"still.unmapped.pre"
grep 'SA:Z:' oot.sam > All_done/$CONTIGS"vs"$REF"odd.mapped.pre" 
grep 'XA:Z:' oot.sam >> All_done/$CONTIGS"vs"$REF"odd.mapped.pre"

ls *.pre >test_prefna.in
while read PREFAS
do
cut -f1,10 $PREFAS |sed 's|^|>|'| sed 's|\t|\n|' >$CONTIGS$PREFAS".fna"
done < test_prefna.in

awk '$2 != 4' oot.sam |awk '$2 <= 200' | cut -f1,4,6| sort -g -k2 > ootcut.sam

cat -v scaffs/$CONTIGS* >test.in

cut -f1 ootcut.sam > testscaffs.txt
	while read SCAFF
	do
echo $SCAFF
echo $SCAFF > test.name
grep -w $SCAFF ootcut.sam | cut -f2,3 > test.stuff
grep -w $SCAFF test.in | cut -f2 | wc -m > test.len

paste test.name test.stuff test.len >test.oot
cat testAll.oot test.oot > temp.arse
mv temp.arse testAll.oot

	done <testscaffs.txt
	
mv testAll.oot ordered/$CONTIGS"orderedto"$REF".txt"

done < contigs.txt


rm -f ref*
rm -f *.arse

#you need to remove (non soft-clipped) scaffold that map within other scaffolds.