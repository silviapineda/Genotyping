
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


setwd("/Users/Pinedasans/Catalyst/Results/Validation/")

load("/Users/Pinedasans/Catalyst/Data/Genotyping/Genotyping_QC.Rdata")
non.list<-seq(1,1326,2)


###To pass a QC on MAF
load("/Users/Pinedasans/Catalyst/Data/Genotyping/MAF_HEW.Rdata")
p.value_hwe_adj<-p.adjust(p.value_hwe,"bonferroni") #389,637 with fdr 131,558 with bonferroni pass the MT correction and are not in HWE
table(maf<0.01) #12,739 has a maf<0.01 and 138,999 has a maf<0.05

id.maf<-match(names(maf),colnames(SNP_calls_diff2))
SNP_calls_diff2_maf<-SNP_calls_diff2[,na.omit(id.maf)]
maf2<-maf[which(is.na(id.maf)==F)]
SNP_calls_diff_mafQC<-SNP_calls_diff2_maf[,which(maf2>0.01)] #801,977

id.snp<-match(colnames(SNP_calls_diff_mafQC),annot$SNP_id)
annot_snp_mafQC<-annot[na.omit(id.snp),]


############################
#### Validation with RF ###
############################
library("randomForest")
library("RColorBrewer")

###########################
###results from Fisher ####
###########################
resultsFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTest.txt",header=T,sep="\t")

#Variants that overlap
id.snp_rs<-match(resultsFisher$snp138,annot_snp_mafQC$SNP_rs)
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.snp_rs)]

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "19variants Overlap")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

###Genes that Overlap
id.gene<-match(resultsFisher$Gene.refGene,annot_snp_mafQC$gene.name)
id.snp<-match(annot_snp_mafQC$SNP_id[id.gene],colnames(SNP_calls_diff_mafQC))
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.snp)] #108 variants
rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "Genes Overlap Fisher")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

#####Plotting MDS 
d <- dist(SNP_rf_selected) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- fit$points[,1]
y <- fit$points[,2]

##by endpoint
COLOR==brewer.pal(3,"Set1")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",col=COLOR[factor(annot_samplePaired$Outcome[non.list])], 
     main="19 variants",pch=20)
legend("topright", legend=levels(annot_samplePaired$Outcome[non.list]), col=COLOR,pch=20)

##by race mismatch
COLOR==brewer.pal(3,"Set1")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",col=COLOR[factor(annot_samplePaired$race_mismatch[non.list])], 
     main="19 variants",pch=20)
legend("topright", legend=c("race match","race mismatch"), col=COLOR,pch=20)

##by center
COLOR=brewer.pal(21,"Set1")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",col=COLOR[factor(annot_samplePaired$Center[non.list])], 
     main="19 variants",pch=20)
legend("topright", legend=levels(factor(annot_samplePaired$Center[non.list])), col=COLOR,pch=20)

####by batch
library(crlmm)
rownames(SNP_calls_diff_mafQC)
plate <- substr(rownames(SNP_calls_diff_mafQC), 12, 13) #1

COLOR=brewer.pal(21,"Set1")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",col=COLOR[factor(plate)], 
     main="19 variants",pch=20)
legend("topright", legend=levels(factor(plate)), col=COLOR,pch=20)

####by HLA region
#Variants that overlap
id.hla<-grep("HLA",annot_snp_mafQC$gene.name)
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.hla)]

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "19variants Overlap")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

##by race mismatch
d <- dist(SNP_rf_selected) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- fit$points[,1]
y <- fit$points[,2]

COLOR==brewer.pal(3,"Set1")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",col=COLOR[factor(annot_samplePaired$race_mismatch[non.list])], 
     main="19 variants",pch=20)
legend("topright", legend=c("race match","race mismatch"), col=COLOR,pch=20)


#######################
##  results from RF ##
######################
resultsRF<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",header=T,sep="\t")

#Variants that overlap
id.snp_rs<-match(resultsRF$snp138,annot_snp_mafQC$SNP_rs)
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.snp_rs)]

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "11variants Overlap")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)


###Genes Overlap
id.gene<-match(resultsRF$Gene.refGene,annot_snp_mafQC$gene.name)
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.gene)] #56 variants
rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "Genes Overlap RF")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

###Randomly 108 variants 
SNP_rf_selected<-SNP_calls_diff_mafQC[,sample(ncol(SNP_calls_diff_mafQC),108)] #108 variants
rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "108 random variants")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)


################################################
##RF with  the variants that are sign p<0.001 ##
################################################

###Runnning Logistic Regression to find significant variants
OR<-rep(NA,ncol(SNP_calls_diff_mafQC))
p.value.OR<-rep(NA,ncol(SNP_calls_diff_mafQC))
for (i in 1:ncol(SNP_calls_diff_mafQC)){
  print(i)
  tab<-table(SNP_calls_diff2[,i])
  if(length(tab)>1){
    model<-glm(annot_samplePaired$Outcome[non.list]~SNP_calls_diff_mafQC[,i],family = "binomial")
    OR[i]<-exp(coef(model))[2]
    p.value.OR[i]<-coefficients(summary(model))[2,4]
  }
}

write.table(cbind(OR,p.value.OR),"logisticRegressionMAFQC.txt")
results_reg<-read.table("logisticRegressionMAFQC.txt")
rownames(results_reg)<-colnames(SNP_calls_diff_mafQC)
results_reg_sign<-results_reg[results_reg$p.value.OR<0.001,] ##1174
id.snp<-match(rownames(results_reg_sign),colnames(SNP_calls_diff_mafQC))
SNP_rf_selected<-SNP_calls_diff_mafQC[,na.omit(id.snp)]

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "SignificantVariants")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)


#############
### Apply this variant in the ExomeSeq data ###
###############
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq_Diff_demo.Rdata")
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq/ExomeSeqVCF_SNPs.Rdata")
demographics<-read.table("/Users/Pinedasans/Catalyst/Data/Demographics.txt",sep="\t",header=T)
non.list<- seq(1,56,2) ##Donors
exome_variants_diff_complete<-exome_variants_diff2[ , ! apply( exome_variants_diff2, 2 , function(x) any(is.na(x)) ) ] #450,981

id.snp<-match(rownames(results_reg_sign),annot_snp_mafQC$SNP_id)
id.snp_exome<-match(annot_snp_mafQC[id.snp,2],df_joint_qc$snp138)
id.snp_id<-match(df_joint_qc[na.omit(id.snp_exome),14],colnames(exome_variants_diff_complete))
exome_variants_rf_selected<-exome_variants_diff_complete[,na.omit(id.snp_id)] 
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected),proximity=TRUE, keep.forest=T,ntree=100)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "Original")
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20)








##Count number of variants per pair that are DIFFERENT
variant_mismatch <- NULL
for (i in 1:511) {
  print(i)
  variant_mismatch[i]<- table(SNP_calls_diff_mafQC[i,]==1)[2]
}
write.table(variant_mismatch,"/Users/Pinedasans/Catalyst/Results/Validation/variant_mismatch_MAFfilter.txt")

boxplot(variant_mismatch~annot_samplePaired$Outcome[non.list],frame.plot = FALSE,col=c("goldenrod","darkolivegreen4"),
        ylab="Num Variants Differ")
anova(lm(variant_mismatch~annot_samplePaired$Outcome[non.list])) #p-value = 0.13

##ICAM1-rs5030351-SNP_A-8348664
variantMismatch<-SNP_calls_diff2[,509897]

samples_mismatch<-annot_samplePaired[match(names(xx[which(xx==1)]),annot_samplePaired$CEL.file),]
SNP_mismatch<-SNP_calls_paired[509897,match(samples_mismatch$CEL.file,colnames(SNP_calls_paired))]
table(samples_mismatch$Outcome,as.numeric(SNP_mismatch)) # 1 and 2 means that the recipient has the variant and 3 that the donor has the variant

outcome<-ifelse(annot_samplePaired$Outcome[non.list]=="TX",0,1)
model<-glm(outcome~variantMismatch,family = "binomial")
exp(coef(model))[2]
coefficients(summary(model))[2,4]
library(effects)
plot(allEffects(model),main="Effect Size VariantMismatch and Outcome")

##AGGF1-rs73125845
id.gene<-grep("AGGF1",annot$gene.name)
id.snp<-match(annot[id.gene,1],colnames(SNP_calls_diff2))
SNP_AGGF1<-SNP_calls_diff2[,na.omit(id.snp)]
OR<-rep(NA,ncol(SNP_AGGF1))
p.value.OR<-rep(NA,ncol(SNP_AGGF1))
for (i in 1:ncol(SNP_AGGF1)){
  print(i)
  tab<-table(SNP_AGGF1[,i])
  if(length(tab)>1){
    model<-glm(annot_samplePaired$Outcome[non.list]~SNP_AGGF1[,i],family = "binomial")
    OR[i]<-exp(coef(model))[2]
    p.value.OR[i]<-coefficients(summary(model))[2,4]
  }
}


###To look for complete cases to run PCA and RF
SNP_calls_diff_complete<-SNP_calls_diff2[ , ! apply( SNP_calls_diff2, 2 , function(x) any(is.na(x)) ) ] ##only 275,869 are complete cases
save(SNP_calls_diff_complete,annot_samplePaired,file="Validation_RF.Rdata")

















##Read variants from the discovery set
resultsFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTest.txt",header=T,sep="\t")

results_reg<-read.table("logisticRegressionMAFQC.txt")
rownames(results_reg)<-colnames(SNP_calls_diff_mafQC)
id.annot<-match(rownames(results_reg),annot$SNP_id)
annot_snps<-annot[na.omit(id.annot),]

##SNps that overlap between both datasets
id.snp<-match(resultsFisher$snp138,annot_snps$SNP_rs) ##19 SNPs overlap 
results_reg[match(annot_snps[na.omit(id.snp),1],rownames(results_reg)),] ##None Validate

##Genes that overlap
id.gene<-match(resultsFisher$Gene.refGene,annot_snps$gene.name) #108 variants within the genes that overlap
results_reg[match(annot_snps[na.omit(id.gene),1],rownames(results_reg)),]

results_reg_sign<-results_reg[match(annot_sign_results$SNP_id,colnames(SNP_calls_diff2)),]
results_RejvsNorej<-cbind(annot_sign_results,results_reg_sign)


SNP_sign<-colnames(SNP_calls_diff2)[which(results_reg$p.value<0.05)]
id.snp.sign<-match(SNP_sign,annot$SNP_id)
annot_sign<-annot[id.snp.sign,]
id.gene<-match(resultsRF$Gene.refGene,annot_sign$gene.name)
annot_sign_results<-annot_sign[na.omit(id.gene),]
results_reg_sign<-results_reg[match(annot_sign_results$SNP_id,colnames(SNP_calls_diff2)),]
results_RejvsNorej<-cbind(annot_sign_results,results_reg_sign)
write.table(results_RejvsNorej,"ResultsLogisti_RejvsNoRej.txt")


###Variants that overlap results:
id.snp<-match(resultsRF$snp138,annot$SNP_rs) ##No SNP overlaping that validate (3 overlap)
annot_overlap<-annot[na.omit(id.snp),]
id.snp.overlap<-match(annot_overlap$SNP_id,rownames(SNP_calls_paired))
results_reg[id.snp.overlap,]


####################################
####Heatmap with the results#######
##################################
library("pheatmap")
library("RColorBrewer")
id.sign.results<-match(results_RejvsNorej$SNP_id,colnames(SNP_calls_diff2))
SNP_calls_sign<-SNP_calls_diff2[,id.sign.results]

annotation_row = data.frame(
  phenotype = annot_samplePaired$Outcome[non.list])
rownames(annotation_row)<-annot_samplePaired$PairID[non.list]
rownames(SNP_calls_sign)<-annot_samplePaired$PairID[non.list]
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
ann_colors = list (phenotype = c(AR = "goldenrod", TX = "darkolivegreen4"))
pheatmap(na.omit(SNP_calls_sign),cluster_rows = T,color = colorRampPalette(c("white", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors)






###Find the variant in genes that overlap only with the RejvsNorej
id.snp<-match(results_RejvsNorej$SNP_id,colnames(SNP_calls_diff2))
SNP_rf_selected<-SNP_calls_diff2[,na.omit(id.snp)]

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected),proximity=TRUE,na.action=na.roughfix)

COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "Genes Overlap RejvsRej")
legend("bottomright", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

###Selection from the variants with p-value<0.001
SNP_rf_selected<-SNP_calls_diff2[,which(results_reg$p.value<0.001)] #935
SNP_rf_selected_completed<-SNP_rf_selected[ , ! apply( SNP_rf_selected, 2 , function(x) any(is.na(x)) ) ] #342
dataset<-data.frame(SNP_rf_selected_completed)
library("VSURF")
set.seed(2)
fit<-VSURF(x = dataset, y=annot_samplePaired$Outcome[non.list],parallel = TRUE,ncores=30)

rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected_completed),proximity=TRUE,na.action=na.roughfix)
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, annot_samplePaired$Outcome[non.list],palette = COLOR,main = "Genes p<0.001 RejvsNoRej")
legend("bottomleft", legend=levels(annot_samplePaired$Outcome[non.list]),col=COLOR, pch = 20)

