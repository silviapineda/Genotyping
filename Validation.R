
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
SNP_calls_paired[509897,]
xx<-abs(diff(as.numeric(SNP_calls_paired[509897,])))[non.list]
yy<-replace(xx,xx==2,1)
table(yy,annot_samplePaired$Outcome[non.list])


##Read variants from the discovery set
results<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",header=T)

id.overlap<-match(results$snp138,annot$SNP_rs)
annot_overlap<-annot[na.omit(na.omit(id.overlap)),] #14

id.overlap_calls<-match(annot_overlap$SNP_id,rownames(SNP_calls_paired))
SNP_calls_paired_overlap<-SNP_calls_paired[na.omit(id.overlap_calls),] #13

SNP_calls_diff<-apply(SNP_calls_paired,1,function(x) abs(diff(x))[non.list])
SNP_calls_diff2<-apply(SNP_calls_diff,2,function(x) replace(x,x==2,1))

SNP_calls_diff_complete<-SNP_calls_diff2[ , ! apply( SNP_calls_diff2, 2 , function(x) any(is.na(x)) ) ] ##only 275,869 are complete cases
save(SNP_calls_diff_complete,annot_samplePaired,file="Validation_RF.Rdata")

p.value<-rep(NA,ncol(SNP_calls_diff2))
for (i in 1:ncol(SNP_calls_diff2)){
  print(i)
  tab<-table(SNP_calls_diff2[,i],annot_samplePaired$Outcome[non.list])
  if(dim(tab)[1]>1){
    p.value[i]<-fisher.test(tab)$p.value
  }
}

write.table(p.value,"p.value.val.txt")

p.value.adj<-p.adjust(p.value,method="fdr") ###None pass MT correction
SNP_sign<-colnames(SNP_calls_diff2)[which(p.value<0.01)]
id.snp.sign<-match(SNP_sign,annot$SNP_id)
annot_sign<-annot[id.snp.sign,]
id.gene<-match(results$Gene.refGene,annot_sign$gene.name)
annot_sign_results<-annot_sign[na.omit(id.gene),]
p.value[match(annot_sign_results$SNP_id,colnames(SNP_calls_diff2))]
