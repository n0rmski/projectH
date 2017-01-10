###this makes a consensus from homozygous vcf (don't use for hets!)###

args <- commandArgs(TRUE)
vcf_file <- read.table(args[1], quote="\"") #arg 1 input file

genotype<-sapply(strsplit(as.character(vcf_file$V10), ":"), "[", 1)
depth<-sapply(strsplit(as.character(vcf_file$V8), ";"), "[", 1)
depth<-as.numeric(gsub("DP=","",depth))
genesample<-data.frame(cbind(as.character(vcf_file$V2),as.character(vcf_file$V4),as.character(vcf_file$V5),genotype, as.character(depth)),stringsAsFactors=FALSE)
row.names(genesample)<-genesample$V1
names(genesample)<-c("position","ref","var","genotype","depth")

genesample$snp1<-sapply(strsplit(as.character(genesample$genotype), "/"), "[", 1)
genesample$snp2<-sapply(strsplit(as.character(genesample$genotype), "/"), "[", 2)
genesample[ is.na(genesample) ] <- 0
genesample$var_called<-as.numeric(genesample$snp1)+as.numeric(genesample$snp2)
genesample$var_called<-ifelse(genesample$var_called>0,1,0)
genesample$var_called<-ifelse(genesample$var_called>0,1,0)

genesample$genotype<-NULL
genesample<-genesample[!genesample$snp1>1, ]
genesample<-genesample[!genesample$snp2>1, ]
genesample$snp1call<-ifelse(genesample$snp1=="0",as.character(genesample$ref), as.character(genesample$var))
genesample$snp2call<-ifelse(genesample$snp2=="0",as.character(genesample$ref), as.character(genesample$var))

geneseq.df<-data.frame(cbind(as.character(genesample$position), as.character(genesample$snp1call)))
names(geneseq.df)<-c("position","call")
write.table(geneseq.df, file=args[2], row.names=F, col.names=T, quote=F)
