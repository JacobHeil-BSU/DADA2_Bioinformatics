## 16S analysis module for Bioinformatics class Day 1
#adapted from the 16S tutorial by Jonas Frankel-Bricker 

##Setup
#required packages
#if you do not have any of these, use the install.packages() function to obtain them
library(dada2)
library(ggplot2)

#set the working directory as the file where this script is contained
setwd("")

#Set a random seed number for the session. It is important to save this number for reproducibility
set.seed(111)

#Set path for raw fastq.gz file location
raw_path <- file.path("fastq_artemisia/")

# Set paths where filtered reads will be deposited. This will generate two new file folders during the filter and trim step and populate them with trimmed reads.
filtpath_R1 <- file.path("fastq_artemisia/filtered_reads/R1/")
filtpath_R2 <- file.path("fastq_artemisia/filtered_reads/R2/")

#Create an object that lists the names of the forward (R1) and reverse (R2) fastq.gz files
fnFs_16S <- sort(list.files(raw_path, pattern = "R1", full.names = TRUE))
fnRs_16S <- sort(list.files(raw_path, pattern = "R2", full.names = TRUE))

#Designate locations of filtered reads
filts_F <- list.files(filtpath_R1, pattern = "R1", full.names = TRUE)
filts_R <- list.files(filtpath_R2, pattern = "R2", full.names = TRUE)

#Consolidate sample names
sample.names_F <- sapply(strsplit(basename(filts_F), "_"), `[`, 2)
sample.names_R <- sapply(strsplit(basename(filts_R), "_"), `[`, 2)
names(filts_F) <- sample.names_F
names(filts_R) <- sample.names_R
dds <- vector("list", length(sample.names_F))
names(dds) <- sample.names_F

##Filter and Trim the raw reads
#First determine the length to which you want to trim the reads
#Next plot the quality profile of all forward and reverse reads
#You can view these plots full screen by clicking the zoom button above the plot in the GUI
plotQualityProfile() + scale_x_continuous(breaks = seq(from = 0, to = 300, by = 25)) + geom_hline(yintercept = 20, color = "red", linetype = 2)

plotQualityProfile() + scale_x_continuous(breaks = seq(from = 0, to = 300, by = 25)) + geom_hline(yintercept = 20, color = "red", linetype = 2)

#input length to which you want to trim both forward and reverse as the parameter trunClen(forward trim length, reverse trim length)
filterAndTrim(fwd = fnFs_16S, filt = , 
              rev = fnRs_16S, filt.rev = , 
              truncLen = c( , ), compress = TRUE, verbose = TRUE)

#Learn errors
errF <- learnErrors(filtpath_R1, multithread = TRUE, verbose = TRUE)
errR <- learnErrors(filtpath_R2, multithread = TRUE, verbose = TRUE)

# for loop to dereplicate sequences, account for estimated errors, merge reads, and produce ASVs
for(sam in sample.names_F){ #the loop repeats once for each name in sample.names_F
  cat("Processing:", sam, "\n") #prints the name of the sample that is currently being run
  derepFs <- derepFastq(filts_F[[sam]]) #dereplicate the forward read for this sample
  derepRs <- derepFastq(filts_R[[sam]]) #dereplicate the reverse read for this sample
  ddf <- dada(derepFs, err = errF) #produce ASVs from forward reads
  ddr <- dada(derepRs, err = errR) #produce ASVs from reverse reads
  merged <- mergePairs(ddf, derepFs, ddr, derepRs, minOverlap = 12) #merge forward and reverse reads
  dds[[sam]] <- merged #store final product in an object with other samples
}

#Create sequence table and remove chimeric sequences
seqtab <- makeSequenceTable()
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus")
View(seqtab_nochim)
