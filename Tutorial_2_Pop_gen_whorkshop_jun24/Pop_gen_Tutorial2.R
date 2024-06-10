###########################################################
# Population Genomic analysis tutorial from the workshop held at Kasetsart university on 10th June 2024. 
# For details on files, see READ me page. We will be working in R studio for the entire workshop.

#Please open R studio and clean your area (optional)

rm(list = ls())

# Create a folder directory for all your files and set as working directory:
setwd("~/Desktop/Tutorial_2_Pop_gen_whorkshop_jun24")

#Upload libraries
library(tidyverse)
library(ggplot2)


#### PCA Ptychadena erlangueri

## Import PCA matrix
pca <- read_table("./Perlangeri_PCA.eigenvec", col_names = FALSE)
pca
## the output of Plink 1.9 for PCA analysis return 2 coluns for samples names
# Let's exclude one
pca <- pca[,-1]
pca
# Rename first column to contain Individuals
names(pca)[1] <- "ind"

# Rename the other columns with PC and number
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca

# import data on percentage of variation
pcaeigenval <- scan("./Perlangeri_PCA.eigenval")

# Plot variance explained
pve <- data.frame(PC = 1:13, pve = pcaeigenval/sum(pcaeigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
pve

# Now let's plot the PCA analysis in 2D

# Plot without grouping info

# Plot without grouping info

PCA_plot1 <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3)
PCA_plot1
#ad information on variation explaine in graph
PCA_plot1 + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


# Plot East and West samples

pop=c("West", "East", "East", "East", "East", "East", "West", "West", "West", "West", "West", "West", "West")

PCA_plot2 <- ggplot(pca, aes(PC1, PC2, label = TRUE, color=pop)) + geom_point(size = 3)
PCA_plot2
PCA_plot2 + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# Plot according mtDNA clades

pop=c("E3-West", "E1-Kibre", "E1-Kibre", "E2-Assela", "E2-Assela", "E2-Assela", "E3-West", "E3-West", "E3-West", "E3-West", "E3-West", "E3-West", "E3-West")

PCA_plot <- ggplot(pca, aes(PC1, PC2, label = TRUE, color=pop)) + geom_point(size = 3)
PCA_plot
PCA_plot + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


#### Fst Ptychadena erlangueri

#Fst E1-E2
Fst_E1_E2<-read.table("popE1_vs_E2_10k_w_FST.windowed.weir.fst", sep="\t", header=T)
head(Fst_E1_E2)
#Removed rows containing missing values
Fst_E1_E2_noNA <- na.omit(Fst_E1_E2)
#summary and plot
summary(Fst_E1_E2_noNA$WEIGHTED_FST)
hist(Fst_E1_E2_noNA$WEIGHTED_FST, col ='azure3')

#Fst E1-E3
Fst_E1_E3<-read.table("popE1_vs_E3_10k_w_FST.windowed.weir.fst", sep="\t", header=T)
head(Fst_E1_E3)
#Removed rows containing missing values
Fst_E1_E3_noNA <- na.omit(Fst_E1_E3)
#summary and plot
summary(Fst_E1_E3_noNA$WEIGHTED_FST)
hist(Fst_E1_E3_noNA$WEIGHTED_FST, col ='green')

#Fst E2-E3
Fst_E2_E3<-read.table("popE2_vs_E3_10k_w_FST.windowed.weir.fst", sep="\t", header=T)
head(Fst_E2_E3)
#Removed rows containing missing values
Fst_E2_E3_noNA <- na.omit(Fst_E2_E3)
#summary and plot
summary(Fst_E2_E3_noNA$WEIGHTED_FST)
hist(Fst_E2_E3_noNA$WEIGHTED_FST, col ='blue')


### Nucleotyde diversity Ptychadena erlangueri populations

#Pi_all
Pi_all<-read.table("Pi_all_pops.windowed.pi", sep="\t", header=T)
head(Pi_all)
#summary and plot
summary(Pi_all$PI)
hist(Pi_all$PI, col ='azure3')

#Pi E1
Pi_E1 <-read.table("Pi_E1.windowed.pi", sep="\t", header=T)
#summary and plot
summary(Pi_E1$PI)
hist(Pi_E1$PI, col ='green')

#Pi E2
Pi_E2 <-read.table("Pi_E2.windowed.pi", sep="\t", header=T)
#summary and plot
summary(Pi_E2$PI)
hist(Pi_E2$PI, col ='blue')

#Pi E3
Pi_E3 <-read.table("Pi_E3.windowed.pi", sep="\t", header=T)
#summary and plot
summary(Pi_E3$PI)
hist(Pi_E3$PI, col ='coral')


###Tajima's D

### E1
TajimaD_E1<-read.table("tajD_10k_E1.Tajima.D", sep="\t", header=T)
#summary and plot
summary(TajimaD_E1$TajimaD)
hist(TajimaD_E1$TajimaD, col ='azure3')

### E2
TajimaD_E2<-read.table("tajD_10k_E2.Tajima.D", sep="\t", header=T)
#summary and plot
summary(TajimaD_E2$TajimaD)
hist(TajimaD_E2$TajimaD, col ='green')

### E3
TajimaD_E3<-read.table("tajD_10k_E3.Tajima.D", sep="\t", header=T)
#summary and plot
summary(TajimaD_E3$TajimaD)
hist(TajimaD_E3$TajimaD, col ='blue')
