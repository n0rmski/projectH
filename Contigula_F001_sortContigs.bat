
ls contigsHERE/*.fna | sed 's|contigsHERE/||g' | sed 's|.fna||' >contigs.txt

rm -f *.fna
rm -f *.pre
while read CONTIGS
do
sed "s|tom|$CONTIGS|" contigsHERE/$CONTIGS".fna" > test.fna

k=101
bwa mem -k$k /home/paul/runfolder/filters/hg19/hg19.fa test.fna -M >test.sam
#bwa mem -k$k /home/paul/runfolder/filters/hg19plusAlt/hg19plusAlt.fa contigs/$CONTIGS".fna" -M >test.sam

java -Xmx2g -jar /home/paul/zprogz/picard-tools-1.129/picard.jar SamFormatConverter INPUT=test.sam OUTPUT=test.bam
samtools sort test.bam test.sorted
java -Xmx2g -jar /home/paul/zprogz/picard-tools-1.129/picard.jar SamFormatConverter INPUT=test.sorted.bam OUTPUT=test.sam

grep 'SA:Z:' test.sam |awk '$2 <= 200' >test.doubledippers.pre
grep 'XA:Z:' test.sam >>test.doubledippers.pre

grep -v "@SQ" test.sam|grep -v "@PG"|awk '$2 <= 200'|grep -v 'XA:Z:' >test2.sam

grep 'XA:Z:chr6' test.sam >>test2.sam

awk '$2 == 4' test2.sam > testkeep.unmapped.pre
awk '$2 != 4' test2.sam > test.mapped.sam

awk '$3 != "chr6"' test.mapped.sam >test.notHLA.pre
awk '$3 == "chr6"' test.mapped.sam >test.chr6.sam

awk '$4 <= 28477896' test.chr6.sam > test.chr6notHLA.pre
awk '$4 >= 33500000' test.chr6.sam >> test.chr6notHLA.pre
awk '$4 >= 28477896' test.chr6.sam | awk '$4 <= 33500000' > test.hla.sam

awk '$5 != 60' test.hla.sam > testkeep.maynotbehla.pre
awk '$5 == 60' test.hla.sam > testkeep.hla.pre

#cat *.pre > testIhopethisisithesameasinput.pre

cat testkeep* > testgroups1to3.pre

ls *.pre >test_prefna.in
while read PREFAS
do
cut -f1,10 $PREFAS |sed 's|^|>|'| sed 's|\t|\n|' >$CONTIGS$PREFAS".fna"
done < test_prefna.in

cut -f1,10 testgroups1to3.pre >forLG/$CONTIGS"groups1to3_forLG.fna"

cp $CONTIGS"testgroups1to3.pre.fna" contigs/
mv $CONTIGS* All_pre/

#rm -f test*
done < contigs.txt
