
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
non.list<-seq(1,1022,2)

load("/Users/Pinedasans/Catalyst/Data/Genotyping/MAF_HEW.Rdata")
p.value_hwe_adj<-p.adjust(p.value_hwe,"bonferroni") #389,637 with fdr 131,558 with bonferroni pass the MT correction and are not in HWE
table(maf<0.05) #12,739 has a maf<0.01 and 138,999 has a maf<0.05

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

###To look for complete cases to run PCA and RF
SNP_calls_diff_complete<-SNP_calls_diff2[ , ! apply( SNP_calls_diff2, 2 , function(x) any(is.na(x)) ) ] ##only 275,869 are complete cases
save(SNP_calls_diff_complete,annot_samplePaired,file="Validation_RF.Rdata")


##Read variants from the discovery set
resultsRF<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF_RejvsNoRej.txt",header=T)

OR<-rep(NA,ncol(SNP_calls_diff2))
p.value.OR<-rep(NA,ncol(SNP_calls_diff2))
for (i in 1:ncol(SNP_calls_diff2)){
  print(i)
  tab<-table(SNP_calls_diff2[,i])
  if(length(tab)>1){
    model<-glm(endpoint~SNP_calls_diff2[,i],family = "binomial")
    OR[i]<-exp(coef(model))[2]
    p.value.OR[i]<-coefficients(summary(model))[2,4]
  }
}

write.table(cbind(OR,p.value),"logisticRegression.txt")

results_reg<-read.table("logisticRegression.txt")
rownames(results_reg)<-colnames(SNP_calls_diff2)

p.value.adj<-p.adjust(results_reg$p.value,method="fdr") ###None pass MT correction
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



###########################################
#### Apply RF in the Sign<0.05 dataset ###
#########################################
library("randomForest")
SNP_rf_selected<-SNP_calls_diff2[,which(results_reg$p.value<0.05)]
SNP_rf_selected_completed<-SNP_rf_selected[ , ! apply( SNP_rf_selected, 2 , function(x) any(is.na(x)) ) ] #14,934
rf_output_total <- randomForest(annot_samplePaired$Outcome[non.list]~.,data=data.frame(SNP_rf_selected_completed),proximity=TRUE)

