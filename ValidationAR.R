
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
library("randomForest")
library("RColorBrewer")


setwd("/Users/Pinedasans/Catalyst/Results/Validation/")

load("/Users/Pinedasans/Catalyst/Data/Genotyping/SNP_calls_diff_imputed.Rdata")
SNP_calls_diff_imputed<-SNP_calls_diff_imputed[,-1]
load("/Users/Pinedasans/Catalyst/Data/Genotyping/SNP_calls_diff.Rdata")
non.list<-seq(1,1326,2)

###Take only the AR
annot_samplePaired_AR<-annot_samplePaired$Outcome[non.list][which(annot_samplePaired$Outcome[non.list]=="AR")]
SNP_calls_diff_AR<-SNP_calls_diff_imputed[which(annot_samplePaired$Outcome[non.list]=="AR"),]

###########################
###results from Fisher ####
###########################
resultsFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTest.txt",header=T,sep="\t")

#############################
##  Variants that overlap ##
###########################
merge_results<-merge(resultsFisher,annot_snp_mafQC,by.x = c("Chr","Start"),by.y = c("Chr","Pos"))
id.snp<-match(merge_results$SNP_id,colnames(SNP_calls_diff_AR))
SNP_rf_selected<-SNP_calls_diff_AR[,id.snp]

#####################################################
### Apply clustering to find the two groups in AR ##
###################################################
d <- dist(as.matrix(SNP_rf_selected))   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc)     

grp <- cutree(hc, k = 2)
COLOR=brewer.pal(3,"Set1")
plot(hc, cex = 0.2) # plot tree
rect.hclust(hc, k = 2, border = COLOR) # add rectangle


############################################################
##  Count number of variants per pair that are DIFFERENT  ##
############################################################
id.AR<-match(names(grp),annot_samplePaired$CEL.file)
annot_samplePaired$Outcome2<-as.character(annot_samplePaired$Outcome)
annot_samplePaired$Outcome2[(id.AR-1)]<-grp
annot_samplePaired$Outcome2[id.AR]<-grp
annot_samplePaired$Outcome2<-replace(annot_samplePaired$Outcome2,annot_samplePaired$Outcome2==1,"AR-1")
annot_samplePaired$Outcome2<-replace(annot_samplePaired$Outcome2,annot_samplePaired$Outcome2==2,"AR-2")

variant_mismatch <- NULL
for (i in 1:663) {
  print(i)
  variant_mismatch[i]<- table(SNP_calls_diff_mafQC[i,]==1)[2]
}

COLOR=brewer.pal(4,"Set1")
boxplot(variant_mismatch~annot_samplePaired$Outcome2[non.list],frame.plot = FALSE,col=COLOR,
        ylab="Num Variants Differ")
summary(lm(variant_mismatch~relevel(factor(annot_samplePaired$Outcome2[non.list]),ref=1)))
summary(lm(variant_mismatch~relevel(factor(annot_samplePaired$Outcome2[non.list]),ref=2)))


####################################
####Heatmap with the results#######
##################################
library("pheatmap")
annot_samplePaired_AR<-annot_samplePaired[which(annot_samplePaired$Outcome2=="AR-1" |
                                                           annot_samplePaired$Outcome2=="AR-2"),]
non.list<-seq(1,280,2)
annotation_row = data.frame(
  phenotype = annot_samplePaired_AR$Outcome2[non.list])
rownames(annotation_row)<-annot_samplePaired_AR$PairID[non.list]
rownames(SNP_rf_selected)<-annot_samplePaired_AR$PairID[non.list]
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
COLOR=brewer.pal(3,"Set1")
ann_colors = list (phenotype = c("AR-1" = COLOR[1], "AR-2" = COLOR[2]))
pheatmap(na.omit(SNP_rf_selected),cluster_rows = T,cluster_cols = T,color = colorRampPalette(c("white", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors,border_color=F)




