rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Affy 6.0 SNP genotypes Donor/Recipient pairs for validation
###
### CITATION: 
###
### PROCESS: Quality of genotyping
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

library(devtools)
#install_github("rscharpf/crlmmCompendium")
library(crlmmCompendium)

setwd("/Users/Pinedasans/Catalyst/Data/Genotyping/")
load("clrm.Rdata")
load("cnSet.Rdata")

load("Genotyping_QC.Rdata")
samples_split<-strsplit(sampleNames(cnSet), "[.]")
samples_names<-unlist(lapply(samples_split, `[[`, 1))
SampleSet <- cnSet[,match(colnames(SNP_calls_paired), samples_names)]

invisible(open(SampleSet))
genotypeSet <- SampleSet[match("SNP_A-8348664", featureNames(SampleSet)), ]
invisible(close(SampleSet))

#extracts the normalized intensities, the genotype calls, and the confidence scores for the genotypes 
#and stores the results in an object of class ‘data.frame’
df <- prePredictPanel(genotypeSet)
fill1 <- brewer.pal(3, "Set1")[df$gt]


##grey scale for the confidence score
gt.conf <- df$gt.conf
min.conf <- min(gt.conf)
max.conf <- max(gt.conf)
sc<-gt.conf>0.9
#sc <- (gt.conf - min.conf)/(max.conf - min.conf)
fill2 <- sapply(sc, grey)

#subset of samples to highlight the batch effects frequently observed in large studies
dt <- strftime(protocolData(genotypeSet)$ScanDate, "%Y-%m-%d",usetz = FALSE)
range(dt)

dt.batch <- split(dt, batch(genotypeSet))
sapply(dt.batch, range)

##select different colors for batch A and batch B, and use white for the remaining samples.
batch.scale <- which(batch(genotypeSet) == "11")
batch.sloth <- which(batch(genotypeSet) == "25")
plate.cols <- brewer.pal(8, "Accent")[c(3, 8)]
fill3 <- rep("white", nrow(df))
fill3[batch.scale] <- plate.cols[1]
fill3[batch.sloth] <- plate.cols[2]

##replicate the ‘data.frame’ object 3 times, attaching a different set of fill colors for each replicate
df2 <- rbind(df, df, df)
df2$fill <- c(fill1, fill2, fill3)
colorby <- c("genotype", "confidence score", "plate")
df2$colorby <- factor(rep(colorby, each = nrow(df)), levels = colorby,ordered=TRUE)

#xyplot creates the desired multi-panel
png("ABfig.png",width = 900, height = 300)
xyplot(A ~ B | colorby, df2, panel = function(x, y,col, fill, plate.cols, ..., subscripts){
  panel.grid(h = 5, v = 5)
  panel.xyplot(x, y, col = "grey60", fill = fill[subscripts], ...)
  if(panel.number() == 3){
    lpoints(x[batch.scale], y[batch.scale], fill = plate.cols[1], ...)
    lpoints(x[batch.sloth], y[batch.sloth], fill = plate.cols[2], ...)
  }
}, aspect = "iso", fill = df2$fill, col = df2$col, cex = 0.6,
pch = 21,plate.cols = plate.cols, xlab = expression(log[2](I[B])),
ylab = expression(log[2](I[A])), main = featureNames(genotypeSet),
layout = c(3,1), par.strip.text = list(lines = 0.9, cex = 0.6))
dev.off()
