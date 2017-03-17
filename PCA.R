rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Esome Sequencing Donor/Recipient pairs
###
### CITATION: 
###
### PROCESS: Principal Components Analysis
###           
### DESCRIP: ExomeSeq Data from Donors/Recipients pairs 
###           27 pairs + 1 pair with two donors
###         
###
### Author: Silvia Pineda
### Date: June, 2016
############################################################################################


#Work Directory
setwd("/Users/Pinedasans/Catalyst/Results/Genotyping/")

###Reading the Total genotypes data frame
load("/Users/Pinedasans/Catalyst/Data/Genotyping/Genotyping_QC.Rdata")

################################
####Principal Component analysis
###############################

##############################################
##   1. Analysis considering only my data ####
##############################################
SNP_calls_diff_complete<-t(SNP_calls_paired[! apply( SNP_calls_paired, 1 , function(x) any(is.na(x)) ), ]) #275,869
write.table()

pca <- prcomp(SNP_calls_diff_complete) 
plot(pca, type = "l")
biplot(pca,cex=c(0.6,0.6))

SPP <- annot_samplePaired$Outcome
#SPP <-addNA(SPP)
#levels.SPP <- factor(c("Rej","NoRej"))
COLOR <- c(2:3)

pc <- c(1,2)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=COLOR[SPP],pch=20,xlab="PCA1",ylab="PCA2")
legend(-90,35, legend=levels(SPP), col=COLOR,pch=20,cex=0.8)



##########################################################
##    2. Analysis considering the 1000G population ######
#########################################################

#########################
####Read the 1000G vcf 
#########################
#Work Directory
setwd("/Users/Pinedasans/Data/Catalyst/ExomeSeq/1000G/")

library("VariantAnnotation") #load the package
variants.1000G.vcf <- readVcf("1000G.variants.selected.vcf","hg19")


library(SNPRelate)
###Convert into a data frame
snpgdsVCF2GDS("1000G.variants.selected.vcf", "1000G.variants.selected.vcf.gds", method="copy.num.of.ref")
genofile <- snpgdsOpen("1000G.variants.selected.vcf.gds")
df.1000G <- data.frame(
  snp_id = rownames(variants.1000G.vcf),
  snp_position = read.gdsn(index.gdsn(genofile, "snp.position")),
  snp_chromosome = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
  snp_allele = read.gdsn(index.gdsn(genofile, "snp.allele")),
  snp_annot = read.gdsn(index.gdsn(genofile,  path="snp.annot/filter"))
)


df.1000G <- cbind(df.1000G,t(read.gdsn(index.gdsn(genofile, "genotype"))))
colnames(df.1000G)[6:dim(df.1000G)[2]]<-colnames(variants.1000G.vcf)
save(df.1000G,file="ExomeSeq1000GVCF.Rdata")

write.table(df.1000G,file="df.1000G.txt")
##I have executed the program Mikel has sent to me (silvia-project-assembly-1.0.jar) to chande the 0 for 2 and 2 for 0.
##The name is df.1000G.mod.txt
df_1000G_mod<-read.table("/Users/Pinedasans/Data/Catalyst/ExomeSeq/1000G/df_1000G_mod.txt") #433,402
load("/Users/Pinedasans/Data/Catalyst/ExomeSeq/ExomeSeqVCF_SNPs.Rdata")
colnames(df_1000G_mod)[2:3]<-c("Start","Chr")
merge_data <- merge(df_joint_qc,df_1000G_mod, by=c("Chr","Start")) #420,569
exome_variants_1000G<-data.matrix(merge_data[,-c(1:15,72:74)])
rownames(exome_variants_1000G)<-merge_data$snp_id.x
exome_variants_1000G_complete<-na.omit(exome_variants_1000G) #402,169

save(exome_variants_1000G_complete,file="/Users/Pinedasans/Data/Catalyst/ExomeSeq/1000G/ExomeSeq1000G_PCA.Rdata")


load("/Users/Pinedasans/Data/Catalyst/ExomeSeq/pca.Rdata")

sample1000g<-read.table("igsr_samples.txt",header=T,sep="\t")
id.1000G<-match(rownames(pca$x),sample1000g$Sample_name)
demographics<-read.table("/Users/Pinedasans/Data/Catalyst/Demographics.txt",sep="\t",header=T) ##Reading the variant list with the phenotype

###Plotting the result with the super population
SPP <- sample1000g$Superpopulation_code[id.1000G]
SPP <-addNA(SPP)
levels.SPP <- factor(c("AFR","AMR","EAS","EUR","SAS","Transplant"))
COLOR <- c(2:6,1)


pc <- c(1,2)
plot(pca$x[,pc[1]][1:2563], pca$x[,pc[2]][1:2563], col=COLOR[SPP],pch=20,xlab="PCA2",ylab="PCA1")
legend(110,100, legend=levels(levels.SPP), col=COLOR,pch=20,cex=0.8)
#text(pca$x[2505:2559,pc[1]], pca$x[2505:2559,pc[2]],labels=rownames(pca$x)[2505:2559],cex=0.8,col=1)
library("zoom")
zm()

####Plot the pca only with our individuals
demographics$phenotype
SPP<-demographics$phenotype[order(demographics$phenotype)]
levels.SPP <- demographics$X[order(demographics$phenotype)]

non.list<-seq(1,56,2)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
COLOR1 <-NULL
set.seed(1)
COLOR1[non.list[1:14]]<- sample(color, 14)
set.seed(1)
COLOR1[non.list[1:14]+1] <- sample(color, 14)
COLOR2 <-NULL
set.seed(3)
COLOR2[non.list[1:7]]<- sample(color, 7)
set.seed(3)
COLOR2[non.list[1:7]+1] <- sample(color, 7)
COLOR3 <-NULL
set.seed(4)
COLOR3[non.list[1:7]]<- sample(color, 7)
set.seed(4)
COLOR3[non.list[1:7]+1] <- sample(color, 7)


demographics$X

###This is taking for the labels in figure 2 in the ExomeSeqDiffAnalysis_v2_112916.R
label= paste("pair",rownames(variant_list),sep="")

pc <- c(1,2)
pca1<-pca$x[,pc[1]][c(2505:2557,2559,2558,2559)]
pca2<-pca$x[,pc[2]][c(2505:2557,2559,2558,2559)]

labelAMR=c("pair1-D","pair1-R","pair4-D","pair4-R","pair7-D","pair7-R","pair8-D","pair8-R","pair11-D","pair11-R",
           "pair15-D","pair15-R","pair16-D","pair16-R","pair18-D","pair18-R","pair19-D","pair19-R","pair20-D","pair20-R",
           "pair21-D","pair21-R","pair24-D","pair24-R","pair26-D","pair26-R","Pair27-D","pair27-R")
par(mfcol=c(1,3))
plot(pca1[order(demographics$phenotype)][1:28], pca2[order(demographics$phenotype)][1:28], 
     col=COLOR1,pch=19,ylim = c(-80,20),xlim=c(-60,160),xlab="PCA2",ylab="PCA1")
legend(-60,20, legend=labelAMR, col=COLOR1,pch=19,cex=0.8)
#text(pca1[order(demographics$phenotype)][1:28]+1, pca2[order(demographics$phenotype)][1:28]+1,labels=levels.SPP[1:28],cex=0.8,col=1)

labelCMR=c("pair2-D","pair2-R","pair3-D","pair3-R","pair9-D","pair9-R","pair13-D","pair13-R","pair17-D","pair17-R","pair22-D","pair22-R",
           "pair25-D","pair25-R")
plot(pca1[order(demographics$phenotype)][29:42], pca2[order(demographics$phenotype)][29:42], 
   col=COLOR2,pch=19,ylim = c(-80,20),xlim=c(-60,160),xlab="PCA2",ylab="PCA1")
legend(-60,20, legend=labelCMR, col=COLOR2,pch=19,cex=0.8)
#text(pca1[order(demographics$phenotype)][29:42]+1, pca2[order(demographics$phenotype)][29:42]+1,labels=levels.SPP[29:42],cex=0.8,col=1)

labelNR=c("pair5-D","pair5-R","pair6-D","pair6-R","pair10-D","pair10-R","pair12-D","pair12-R","pair14-D","pair14-R","pair23-D","pair23-R",
          "pair28-D","pair28-R")
plot(pca1[order(demographics$phenotype)][43:56], pca2[order(demographics$phenotype)][43:56], 
     col=COLOR3,pch=19,ylim = c(-80,20),xlim=c(-60,160),xlab="PCA2",ylab="PCA1")
legend(-60,20, legend=labelNR, col=COLOR3,pch=19,cex=0.8)
#text(pca1[order(demographics$phenotype)][43:56]+1, pca2[order(demographics$phenotype)][43:56]+1,labels=levels.SPP[43:56],cex=0.8,col=1)

##Obtain a mesuare of distance in the plot
PCA1<-pca$x[2505:2559,pc[1]]
PCA1<-PCA1[c(1:53,55,54,55)]
PCA2<-pca$x[2505:2559,pc[2]]
PCA2<-PCA2[c(1:53,55,54,55)]
non.list<-seq(1,56,2)
x1<-PCA1[non.list]
y1<-PCA2[non.list]
x2<-PCA1[non.list+1]
y2<-PCA1[non.list+1]

Distance <- sqrt(((x2-x1)^2)+((y2-y1)^2))
write.table(Distance,"/Users/Pinedasans/Documents/Catalyst/Results/Distance_PCA_1000G.txt")

x1<-cbind(PCA1[non.list],PCA2[non.list])
x2<-cbind(PCA1[non.list+1],PCA1[non.list+1])
library("fields")
rdist(x1, x2)


###Plot with all the populations
SPP <- sample1000g$Population_name
SPP<-
library(RColorBrewer)
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

COLOR <- col_vector
pc <- c(1,2)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=COLOR[SPP],pch=19)
legend(110,100, legend=levels(SPP), col=COLOR,pch=1,cex=0.8)

