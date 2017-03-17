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

###This has been run in the server!!!!
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

load("clrm.Rdata")
load("cnSet.Rdata")

###Calls
snp.index <- which(isSnp(cnSet) & !is.na(chromosome(cnSet)))
calls<-assayData(cnSet)$call[snp.index,] ## 905422 (SNPs)  1855 (Samples)  genotype calls (1 - AA; 2 - AB; 3 - BB) 
##Delete samples with SNR<5 as a low quality sample
SNR <- cnSet$SNR[]
calls_SNRqc<-calls[,which(SNR>=5)] #1,724 samples remind
####Delete SNPs with bad quality score < 0.25
good_SNPs<-featureNames(crlmmResult)[which(featureData(crlmmResult)$SNPQC>0.25)]
id.good_SNPs<-match(good_SNPs,rownames(calls_SNRqc))
calls_SNRqc_SNPsgood<-calls_SNRqc[na.omit(id.good_SNPs),] ##891,125   1,724

###Genotype confidence Score
confScore<-confs(crlmmResult)
id.calling<-match(rownames(calls_SNRqc_SNPsgood),rownames(confScore))
confScore_call<-confScore[na.omit(id.calling),]
id.sample<-match(colnames(calls_SNRqc_SNPsgood),colnames(confScore_call))
confScore_call_sample<-confScore_call[,na.omit(id.sample)]

calls<-calls_SNRqc_SNPsgood
for (i in 1:ncol(calls)){
  print(i)
  calls[,i]<-replace(calls[,i],confScore_call_sample[,i]<=0.9,NA)
}

save(calls,file="callsMissing.Rdata")

###SNP call rate 95%
calls_snp_crate<-calls[rowSums(is.na(calls)==F)>=1637,] #870,905 after call rate QC

###Sample call rate 95%
calls_genotypes<-calls_snp_crate[,colSums(is.na(calls_snp_crate)==F)>=827359] #1,720 after call rate QC
save(calls_genotypes,file="calls_genotypes.Rdata")

load("/Users/Pinedasans/Catalyst/Data/Genotyping/calls_genotypes.Rdata")

###############################
# Calculate the HWE and MAF ###
###############################
library(HardyWeinberg)
p.value_hwe<-rep(NA,nrow(calls_genotypes))
maf<-rep(NA,nrow(calls_genotypes))
for (i in 1:nrow(calls_genotypes)){
  print(i)
  tab<-table(as.matrix(calls_genotypes[i,]))
  if(length(tab)==3){
    names(tab)<-c("AA","AB","BB")
    maf[i]<-maf(as.vector(tab/sum(tab)))
    hwe<-HWAlltests(tab,verbose=TRUE)
    p.value_hwe[i]<-hwe[1,2]
  }
}

save(maf,p.value_hwe,file="MAF_HEW.Rdata")

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
id.snp.annot<-match(rownames(calls_genotypes),annot$SNP_id)
annot2<-annot[na.omit(id.snp.annot),] 

id.snp.aut<-match(annot2$SNP_id,rownames(calls_genotypes))
calls_genotypes_autosome<-calls_genotypes[na.omit(id.snp.aut),] ##836,872 (final set of SNPs)


####################################################
## Prepare data with the annotation from samples ###
####################################################
annot_sample<-read.table("/Users/Pinedasans/Catalyst/Data/Genotyping/annot_samples.txt",header = T,sep="\t")
samples_split<-strsplit(colnames(calls_genotypes_autosome), "[.]")
samples_names<-unlist(lapply(samples_split, `[[`, 1))
id.samples_calls<-match(annot_sample$CEL.file,samples_names)
SNP_calls<-calls_genotypes_autosome[,na.omit(id.samples_calls)] #1,513 samples (
colnames(SNP_calls)<-annot_sample$CEL.file[which(is.na(id.samples_calls)==F)]
id.samples_annot<-match(colnames(SNP_calls),annot_sample$CEL.file)
annot_sample_2<-annot_sample[id.samples_annot,]
annot_sample_3<-annot_sample_2[order(annot_sample_2$PairID,annot_sample_2$D.R),]
write.table(annot_sample_3,"/Users/Pinedasans/Catalyst/Data/Genotyping/annot_samples_to_delete.txt",row.names = F)

##After deleting the CAN samples
annot_sample_ARTX<-read.table("/Users/Pinedasans/Catalyst/Data/Genotyping/annot_samples_ARTX.txt",header=T,sep="\t")
i=1
annot_samplePaired<-NULL
while (i <=nrow(annot_sample_ARTX)) {
  print(i)
  if(annot_sample_ARTX$PairID[i]==annot_sample_ARTX$PairID[i+1]){
    annot_samplePaired<-rbind(annot_samplePaired,annot_sample_ARTX[i:(i+1),])
    i=i+2
  } else {
      i=i+1
  }
}
###1,022 samples paired

id.snp_paired<-match(annot_samplePaired$CEL.file,colnames(SNP_calls))
SNP_calls_paired<-SNP_calls[,id.snp_paired] #836,872   1,022 final set 

###To obtain the variable with the difference
non.list<-seq(1,1022,2)
SNP_calls_diff<-apply(SNP_calls_paired,1,function(x) abs(diff(x))[non.list])
SNP_calls_diff2<-apply(SNP_calls_diff,2,function(x) replace(x,x==2,1))

save(SNP_calls_paired,SNP_calls_diff2,annot_samplePaired,annot,file="/Users/Pinedasans/Catalyst/Data/Genotyping/Genotyping_QC.Rdata")








# ##########################################
# ### Prepare MAP and PED file for PLINK ###
# ##########################################
# 
# # Create map file
# # third column of mapfile (genetic position (cM) is not in annotation,take one that is mostly 0)
# xx<-rep(0,dim(annot2)[1])
# mapfile<-cbind(annot2[,c("Chr","SNP_rs")],xx,annot2[,c("Pos")])
# colnames(mapfile)[4]<-"Pos"
# 
# 
# # Create ped file
# alleleA<-alleleB<-matrix("", nrow=nrow(calls_SNRqc_SNPsgood_aut), ncol=ncol(calls_SNRqc_SNPsgood_aut),
#                          dimnames=dimnames(calls_SNRqc_SNPsgood_aut))
# # take complement of SNPs on reverse strand
# allA<-annot$Allele.A
# allComplA<-c("A","C","G","T")[match(allA,c("T","G","C","A"))]
# forwA<-ifelse(annot$Strand=="+",allA,allComplA)
# allB<-annot$Allele.B
# allComplB<-c("A","C","G","T")[match(allB,c("T","G","C","A"))]
# forwB<-ifelse(annot$Strand=="+",allB,allComplB)
# for (r in 1:ncol(calls)){
#   alleleA[,r]<-ifelse(calls[,r]<3,forwA,forwB)
#   alleleB[,r]<-ifelse(calls[,r]<2,forwA,forwB)
# }
# ped<-matrix("",ncol=nrow(alleleA)*2+6,nrow=ncol(alleleA))
# for (r in 1:ncol(calls)){
#   ped[r,(1:nrow(alleleA))*2+5]<-alleleA[,r]
#   ped[r,(1:nrow(alleleB))*2+6]<-alleleB[,r]
# }
# # Family id
# ped[,1]<-targets$PatientId
# # sample id
# ped[,2]<-make.names(targets$DisplayLabel)
# # Paternal id, maternal id
# ped[,3:4]<-"0"
# # Gender (1=male, 2=female)
# ped[,5]<-ifelse(targets$Gender=="M","1","2")
# # Affection (0=unknown,1=unaffected, 2=affected)
# ped[,6]<-ifelse(targets$NorTum=="N","1","2")
# #
# write.table(mapfile, file="SNPollier.map", quote=FALSE,
#             row.names=FALSE, col.names=FALSE)
# write.table(ped, file="SNPollier.ped", quote=FALSE, row.names=FALSE,
#             col.names=FALSE)
# 
