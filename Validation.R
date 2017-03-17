
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

##ICAM1-rs5030351-SNP_A-8348664
xx<-SNP_calls_diff2[,509897]
samples_mismatch<-annot_samplePaired[match(names(xx[which(xx==1)]),annot_samplePaired$CEL.file),]
SNP_mismatch<-SNP_calls_paired[509897,match(samples_mismatch$CEL.file,colnames(SNP_calls_paired))]
table(samples_mismatch$Outcome,as.numeric(SNP_mismatch)) # 1 and 2 means that the recipient has the variant and 3 that the donor has the variant


endpoint<-ifelse(annot_samplePaired$Outcome[non.list]=="TX",0,1)
model<-glm(endpoint~yy,family = "binomial")
exp(coef(model))[2]
coefficients(summary(model))[2,4]

SNP_calls_diff_complete<-SNP_calls_diff2[ , ! apply( SNP_calls_diff2, 2 , function(x) any(is.na(x)) ) ] ##only 275,869 are complete cases
save(SNP_calls_diff_complete,annot_samplePaired,file="Validation_RF.Rdata")


##Read variants from the discovery set
results<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF_RejvsNoRej.txt",header=T)

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

p.value.adj<-p.adjust(p.value,method="fdr") ###None pass MT correction
SNP_sign<-colnames(SNP_calls_diff2)[which(p.value<0.05)]
id.snp.sign<-match(SNP_sign,annot$SNP_id)
annot_sign<-annot[id.snp.sign,]
id.gene<-match(results$Gene.refGene,annot_sign$gene.name)
id.snp<-match(results$snp138,annot_sign$SNP_rs) ##No SNP overlaping that validate (3 overl)
annot_sign_results<-annot_sign[na.omit(id.gene),]
p.value.sign<-p.value[match(annot_sign_results$SNP_id,colnames(SNP_calls_diff2))]
OR.sign<-OR[match(annot_sign_results$SNP_id,colnames(SNP_calls_diff2))]
cbind(annot_sign_results,p.value.sign,OR.sign)
results_RejvsNorej<-cbind(annot_sign_results,p.value.sign,OR.sign)
write.table(cbind(annot_sign_results,p.value.sign,OR.sign),"ResultsLogisti_RejvsNoRej.txt")
