
ls contigs/*testgroups1to3.pre.fna | sed 's|contigs/||' |sed 's|testgroups1to3.pre.fna||' >samplist.txt

   while read SAMPLE
   do
echo $SAMPLE
#gunzip Seqs/$SAMPLE*

rm -f ref*

bowtie-build -f -o 1 contigs/$SAMPLE"testgroups1to3.pre.fna" ref

bowtie -q -p28 --trim5 1 --trim3 3 -v 1 -e 10 -I 75 -X 1000 ref -1 Seqs/$SAMPLE"_minusHg19v0_1.fastq" -2 Seqs/$SAMPLE"_minusHg19v0_2.fastq" dump --al test_contigs.fastq --un forLG/$SAMPLE"_Contigsnot.fastq"

#mv *.fastq All_done/

##part 2

echo $SAMPLE
rm -f oot.sam
rm -f ref.fas
cp contigs/$SAMPLE"testgroups1to3.pre.fna" ref.fas

bwa index ref.fas
k=100
bwa mem -t16 -k$k ref.fas forLG/$SAMPLE"_Contigsnot_1.fastq" forLG/$SAMPLE"_Contigsnot_2.fastq" -M >oot.sam


awk '$2 == 113' oot.sam >$SAMPLE"linkers_R1.pre"
awk '$2 == 97' oot.sam >>$SAMPLE"linkers_R1.pre"
awk '$2 == 81' oot.sam >>$SAMPLE"linkers_R1.pre"
awk '$2 == 65' oot.sam >>$SAMPLE"linkers_R1.pre"
awk '$2 == 177' oot.sam >$SAMPLE"linkers_R2.pre"
awk '$2 == 145' oot.sam >>$SAMPLE"linkers_R2.pre" 
awk '$2 == 161' oot.sam >>$SAMPLE"linkers_R2.pre"
awk '$2 == 129' oot.sam >>$SAMPLE"linkers_R2.pre"

mv *.pre forLG/

# rm -f Seqs/$SAMPLE*


  done < samplist.txt
# # # 
