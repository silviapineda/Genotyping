rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Affy 6.0 SNP genotypes Donor/Recipient pairs for validation
###
### CITATION: 
###
### PROCESS: Validation analysis
###           
### DESCRIP: genotypes from Donors/Recipients pairs 
###           
###         
###
### Author: Silvia Pineda
### Date: MArch, 2016
############################################################################################

#source("https://bioconductor.org/biocLite.R")
library(crlmm)
library(genomewidesnp6Crlmm)
library(ff)

library("genefilter")
library("IRanges")
library("MASS")
library("VanillaICE")


setwd("/Users/Pinedasans/Catalyst/Data/Genotyping/")

###Read the CEL files
outdir  <- "/Users/Pinedasans/Catalyst/Data/Genotyping/CEL/"
filenames <- list.celfiles(outdir, full.names = TRUE, pattern = ".CEL")
plate <- substr(basename(filenames), 12, 13) #1885

##Call for genotypes
cnSet <- genotype(filenames = filenames, cdfName = "genomewidesnp6", batch = plate)
save(cnSet,file="/Users/Pinedasans/Catalyst/Data/Genotyping/cnSet.Rdata")
crlmmResult <- crlmm(filenames,cdfName = c("GenomeWideSnp6"))
save(crlmmResult,file="/Users/Pinedasans/Catalyst/Data/Genotyping/clrm.Rdata")

load("/Users/Pinedasans/Catalyst/Data/Genotyping/clrm.Rdata")
load("/Users/Pinedasans/Catalyst/Data/Genotyping/cnSet.Rdata")

###Calls
snp.index <- which(isSnp(cnSet) & !is.na(chromosome(cnSet)))
calls<-assayData(cnSet)$call[snp.index,] ## 905422 (SNPs)  1855 (Samples)  genotype calls (1 - AA; 2 - AB; 3 - BB) 

##Delete samples with SNR<4 as a low quality sample
SNR <- cnSet$SNR[]
calls_SNRqc<-calls[,which(SNR>4)] #1801 samples remind

####Delete SNPs with bad quality score < 0.25
good_SNPs<-featureNames(crlmmResult)[which(featureData(crlmmResult)$SNPQC>0.25)]
id.good_SNPs<-match(good_SNPs,rownames(calls_SNRqc))
calls_SNRqc_SNPsgood<-calls_SNRqc[na.omit(id.good_SNPs),] ##891147   1801

###Genotype confidence Score
confScore<-confs(crlmmResult)
id.calling<-match(rownames(calls_SNRqc_SNPsgood),rownames(confScore))
confScore_call<-confScore[na.omit(id.calling),]
calls_SNRqc_SNPsgood[confScore < 0.9] <- NA

###SNP call rate 95%
calls_snp_crate<-calls_SNRqc_SNPsgood[rowSums(is.na(calls_SNRqc_SNPsgood)==F)>= 1786,]

###Sample call rate 95%
calls_genotypes<-calls_snp_crate[,colSums(is.na(calls_snp_crate)==F)>=846589]



################################
## Calculate the HWE and MAF ###
###############################
library(HardyWeinberg)

p.value_hwe<-NULL
maf<-NULL
for (i in 1:nrow(example_con)){
  print(i)
  tab<-table(calls_SNRqc_SNPsgood[i,])
  names(tab)<-c("AA","AB","BB")
  maf[i]<-maf(as.vector(tab/sum(tab)))
  hwe<-HWAlltests(table(example_con[i,]),verbose=TRUE)
  p.value_hwe[i]<-hwe[1,2]
}




#######################################
## Prepare the SNP annotation file ###
######################################
# library(pd.genomewidesnp.6)
# con <- db(pd.genomewidesnp.6)
# annot <- dbGetQuery(con, "select man_fsetid, dbsnp_rs_id, chrom, physical_pos, strand, gene_assoc from featureSet;")
# annot.noNA<-annot[which(is.na(annot$physical_pos)==F),]
# annot.autosome<-annot.noNA[which(annot.noNA$chrom!="X" & annot.noNA$chrom!="Y" & annot.noNA$chrom!="MT"),]
### Separate the gene information with pyhton program
# write.table(annot.autosome,"/Users/Pinedasans/Catalyst/Data/Genotyping/annotation.txt",row.names = F,col.names = F,sep="\t") 

##Read the annotation file
annot<-read.table("/Users/Pinedasans/Catalyst/Data/Genotyping/annotation_genes.txt")
colnames(annot)<-c("SNP_id","SNP_rs","Chr","Pos","strand","gene.func","gene.name")

## convert strand info to +,-,*
annot$strand[is.na(annot$strand)] <- 2
annot$strand <- annot$strand + 1
annot$strand <- c("+", "-", "*")[annot$strand]

### Match annotation file with SNP file
id.snp.annot<-match(rownames(calls_SNRqc_SNPsgood),annot$SNP_id)
annot2<-annot[na.omit(id.snp.annot),] ##855458 (final set of SNPs)

id.snp.aut<-match(annot2$SNP_id,rownames(calls_SNRqc_SNPsgood))
calls_SNRqc_SNPsgood_aut<-calls_SNRqc_SNPsgood[na.omit(id.snp.aut),] ##855458 (final set of SNPs)


####################################################
## Prepare data with the annotation from samples ###
####################################################
annot_sample<-read.table("/Users/Pinedasans/Catalyst/Data/Genotyping/annot_samples_del.txt",header = T,sep="\t")
samples_split<-strsplit(colnames(SNP_calls_2), "[.]")
samples_names<-unlist(lapply(samples_split, `[[`, 1))
id.samples_calls<-match(annot_sample$CEL.file,samples_names)
SNP_calls_paired<-SNP_calls_2[,na.omit(id.samples_calls)] #1330 samples ( )
colnames(SNP_calls_paired)<-annot_sample$CEL.file[which(is.na(id.samples_calls)==F)]

#id<-match(colnames(SNP_calls_paired),annot_sample$CEL.file)
#annot_sample[-id,] ##This samples have been and rerun above code because there is no info on the recipients


##ICAM1-rs5030351
SNP_calls_paired[526245,]
non.list<-seq(1,1330,2)
xx<-abs(diff(SNP_calls_paired[526245,]))[non.list]
yy<-replace(xx,xx==2,1)
table(yy,annot_sample$Outcome[non.list])


SNP_calls_paired_qc<-SNP_calls_paired[,colSums(is.na(SNP_calls_paired) == FALSE) >= 1197] # samples with a callrate < 90% 






##########################################
### Prepare MAP and PED file for PLINK ###
##########################################

# Create map file
# third column of mapfile (genetic position (cM) is not in annotation,take one that is mostly 0)
xx<-rep(0,dim(annot2)[1])
mapfile<-cbind(annot2[,c("Chr","SNP_rs")],xx,annot2[,c("Pos")])
colnames(mapfile)[4]<-"Pos"


# Create ped file
alleleA<-alleleB<-matrix("", nrow=nrow(calls_SNRqc_SNPsgood_aut), ncol=ncol(calls_SNRqc_SNPsgood_aut),
                         dimnames=dimnames(calls_SNRqc_SNPsgood_aut))
# take complement of SNPs on reverse strand
allA<-annot$Allele.A
allComplA<-c("A","C","G","T")[match(allA,c("T","G","C","A"))]
forwA<-ifelse(annot$Strand=="+",allA,allComplA)
allB<-annot$Allele.B
allComplB<-c("A","C","G","T")[match(allB,c("T","G","C","A"))]
forwB<-ifelse(annot$Strand=="+",allB,allComplB)
for (r in 1:ncol(calls)){
  alleleA[,r]<-ifelse(calls[,r]<3,forwA,forwB)
  alleleB[,r]<-ifelse(calls[,r]<2,forwA,forwB)
}
ped<-matrix("",ncol=nrow(alleleA)*2+6,nrow=ncol(alleleA))
for (r in 1:ncol(calls)){
  ped[r,(1:nrow(alleleA))*2+5]<-alleleA[,r]
  ped[r,(1:nrow(alleleB))*2+6]<-alleleB[,r]
}
# Family id
ped[,1]<-targets$PatientId
# sample id
ped[,2]<-make.names(targets$DisplayLabel)
# Paternal id, maternal id
ped[,3:4]<-"0"
# Gender (1=male, 2=female)
ped[,5]<-ifelse(targets$Gender=="M","1","2")
# Affection (0=unknown,1=unaffected, 2=affected)
ped[,6]<-ifelse(targets$NorTum=="N","1","2")
#
write.table(mapfile, file="SNPollier.map", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(ped, file="SNPollier.ped", quote=FALSE, row.names=FALSE,
            col.names=FALSE)

