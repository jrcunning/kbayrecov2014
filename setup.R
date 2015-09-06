# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning

# =================================================================================================
# DATA PREPARATION
# =================================================================================================
# • Load libraries --------------------------------------------------------------------------------
library(lme4); library(MASS); library(reshape2); library(lattice); library(lmerTest) 
library(LMERConvenienceFunctions); library(pbkrtest); library(scales); library(RColorBrewer) 
library(merTools); library(devtools)
## SPIDA package available at http://r-forge.r-project.org/projects/spida/
#system(paste("svn checkout svn://svn.r-forge.r-project.org/svnroot/spida/"))
#devtools::install("spida/pkg")
library(spida)

addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.8),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}
# -------------------------------------------------------------------------------------------------
# • Load data -------------------------------------------------------------------------------------
# Define data calling function
qPCR <- function(files=list(), sym.target=list(), host.target=NULL, 
                 fluor.norm=list(),
                 sym.copy.number=list(), host.copy.number=1,
                 sym.ploidy=1, host.ploidy=1,
                 sym.extract=1, host.extract=1) {
  require(plyr); require(reshape2)
  data <- do.call("rbind", lapply(files, function(file) {
    data.frame(Filename=file, read.csv(file, skip=6, na.strings="Undetermined")[, 1:17])
  }))
  # Remove rows with no Target specified (empty wells)
  data <- data[which(data$Target.Name!=""), ]
  data$Target.Name <- factor(as.character(data$Target.Name))
  # Code unique sample-plate IDs to distinguish samples run on multiple plates
  data$Sample.Run <- interaction(data$Filename, data$Sample.Name, sep="~")
  # Calculate mean Ct values and SDs of each target for each sample run
  ctmeans <- dcast(data, Sample.Run ~ Target.Name, mean, na.rm=F, value.var="C_")
  ctmeans[is.na(ctmeans)] <- NA
  colnames(ctmeans) <- c(colnames(ctmeans)[1], paste(colnames(ctmeans)[-1], "Ct.mean", sep="."))
  ctsds <- dcast(data, Sample.Run ~ Target.Name, sd, na.rm=T, value.var="C_")
  colnames(ctsds) <- c(colnames(ctsds)[1], paste(colnames(ctsds)[-1], "Ct.sd", sep="."))
  # Calculate number of SD's away from mean of host Ct value for each sample run
  mean.hostCt <- mean(data[which(data$Target.Name==host.target), "C_"], na.rm=T)
  sd.hostCt <- sd(data[which(data$Target.Name==host.target), "C_"], na.rm=T)
  data[which(data$Target.Name==host.target), "sd.hostCt"] <- 
    (data[which(data$Target.Name==host.target), "C_"] - mean.hostCt) / sd.hostCt
  sd.hostCt <- aggregate(data$sd.hostCt, by=list(Sample.Run=data$Sample.Run), FUN=mean, na.rm=T)
  colnames(sd.hostCt)[2] <- paste(host.target, "dist.Ct.mean", sep=".")
  # Combine Ct means, SDs, and host Ct distance
  result <- join_all(list(ctmeans, ctsds, sd.hostCt))
  # Split Sample.Run column into Run (filename) and Sample.Name columns
  result <- cbind(colsplit(as.character(result$Sample.Run), pattern="~", names=c("Filename", "Sample.Name")),
                  result[, -1])
  for (sym in sym.target) {
    sym.SH <- 2^((result[, paste(host.target, "Ct.mean", sep=".")] - fluor.norm[[host.target]]) - 
                   (result[, paste(sym, "Ct.mean", sep=".")] - fluor.norm[[sym]]))
    sym.SH <- sym.SH / (sym.copy.number[[sym]] / host.copy.number)  / (sym.ploidy / host.ploidy) / (sym.extract / host.extract)
    sym.SH[is.na(sym.SH)] <- 0
    result[, paste(sym, "SH", sep=".")] <- sym.SH
  }
  return(result)
}

# Get list of plate files to read in
Mcap.plates <- list.files(path="qPCRdata", pattern="csv$", full.names=T)

# Read in data and calculate S/H ratios
Mcap <- qPCR(Mcap.plates, sym.target=list("C", "D"), host.target="Mcap",
             fluor.norm=list(C=0, D=0, Mcap=0),
             sym.copy.number=list(C=1, D=1), host.copy.number=1,
             sym.ploidy=1, host.ploidy=2,
             sym.extract=0.813, host.extract=0.982)
Mcap[which(Mcap$Sample.Name=="NTC"), ]   # Check NTC's
Mcap <- Mcap[-which(Mcap$Sample.Name=="NTC"), ]  # Remove NTC's

# Parse sample names and dates
sample.names <- rbind.fill(lapply(strsplit(as.character(Mcap$Sample.Name), split="_"), 
                                  function(X) data.frame(t(X))))
colnames(sample.names) <- c("sample", "date", "note")
Mcap <- cbind(sample.names, Mcap[, -2])
Mcap$date <- factor(Mcap$date, levels=c("10.24", "11.04", "11.24", "12.16", "01.14", "05.06"))

# Calculate total S/H ratio and D/C ratio
Mcap$tot.SH <- Mcap$C.SH + Mcap$D.SH
Mcap$logDC <- log(Mcap$D.SH / Mcap$C.SH)

# Identify symbiont clades present (C=C only, CD=C > D, DC=D > C, D=D only)
Mcap$syms <- factor(ifelse(Mcap$C.SH > Mcap$D.SH, ifelse(Mcap$D.SH!=0, "CD", "C"), 
                           ifelse(Mcap$D.SH > Mcap$C.SH, ifelse(Mcap$C.SH!=0, "DC", "D"), NA)),
                    levels=c("C", "CD", "DC", "D"))
# Identify dominant symbiont clade
Mcap$dom <- factor(substr(as.character(Mcap$syms), 1, 1))
# Set zeros to NA to facilitate log transformation
Mcap$C.SH[which(Mcap$C.SH==0)] <- NA
Mcap$D.SH[which(Mcap$D.SH==0)] <- NA
Mcap$tot.SH[which(Mcap$tot.SH==0)] <- NA

table(is.na(Mcap$tot.SH))  # 8 SAMPLES HAVE NO DATA!!!!!!!!!!!!!

# Assign visual ID and reef location metadata
Mcap$vis <- factor(ifelse(as.numeric(as.character(Mcap$sample)) %% 2 == 0, "not bleached", "bleached"))
Mcap$reef <- factor(ifelse(as.numeric(as.character(Mcap$sample)) <= 50, "HIMB",
                           ifelse(as.numeric(as.character(Mcap$sample)) >= 101, "44", "25")))

# Replace date values, create date factor, POSIX date, and numeric days
Mcap$fdate <- revalue(Mcap$date, c("10.24"="20141024", "11.04"="20141104", "11.24"="20141124",
                                   "12.16"="20141216", "01.14"="20150114", "05.06"="20150506"))
Mcap$date <- as.Date(Mcap$fdate, format="%Y%m%d")
Mcap$time <- as.POSIXct(Mcap$date)
Mcap$days <- as.numeric(Mcap$date) - min(as.numeric(Mcap$date))
# -------------------------------------------------------------------------------------------------
# • Filter duplicates (samples run multiple times) --------------------------------------------------
filter.dups <- function(data) {
  keep <- data.frame()  # Create empty data frame to receive runs to keep
  # Identify duplicated samples (i.e., samples run multiple times)
  dups <- unique(rbind(data[duplicated(interaction(data$sample, data$date)), ], 
                       data[rev(duplicated(rev(interaction(data$sample, data$date)))), ]))
  # Create list of data frames for each duplicated sample
  dups.list <- split(dups, f=interaction(dups$sample, dups$date))
  # Loop through each duplicated sample and select which run to keep
  for (sample in dups.list) {
    # If any of the runs are a re-extraction, keep the extraction that gave lowest Mcap Ct value
    if ("reex" %in% sample$note) {  
      # label sample by original extraction or reextraction
      sample$note <- ifelse(Vectorize(isTRUE)(sample$note=="reex"), "reex", "original")
      # get mean Mcap Ct's of original extraction and reextracted sample runs
      meanMcap <- aggregate(sample$Mcap.Ct.mean,  
                            by=list(extract=sample$note), 
                            FUN=mean)
      # choose which extraction to keep based on which has lower mean Mcap Ct
      keepextract <- meanMcap$extract[which.min(meanMcap$x)]  
      if (keepextract=="reex") sample <- sample[which(sample$note=="reex"), ]
      if (keepextract=="original") sample <- sample[which(sample$note=="original"), ]
    }
    # If multiple runs of the original or re-extracted sample still exist, keep the run with
    #   the lowest average standard deviation of technical replicates of each target
    sample <- sample[which.min(rowMeans(sample[, c("Mcap.Ct.sd", "C.Ct.sd", "D.Ct.sd")], na.rm=T)), ]
    # Collect runs to keep
    keep <- merge(keep, sample, all=T)
  }
  # Filter out those runs not kept from data
  nondups <- data[setdiff(rownames(data), rownames(dups)), ]
  result <- merge(nondups, keep, all=T)
  return(result)
}

Mcap.f <- filter.dups(Mcap)

# # Analyze effect of filtering on Mcap Ct values
# par(mfrow=c(2,1), mar=c(3,6,3,1))
# boxplot(Mcap$Mcap.Ct.mean, horizontal=T, ylim=c(20,40), main="All runs", 
#         names="Mcap.Ct.mean", show.names=T, las=2)
# boxplot(Mcap.f$Mcap.Ct.mean, horizontal=T, ylim=c(20,40), main="Filtered runs",
#         names="Mcap.Ct.mean", show.names=T, las=2)
# # Analyze effect of filtering on technical replicate sd values
# par(mfrow=c(2,1), mar=c(3,6,3,1))
# boxplot(Mcap[, c("Mcap.Ct.sd", "C.Ct.sd", "D.Ct.sd")], horizontal=T, las=2, ylim=c(0,5),
#         main="All runs")
# boxplot(Mcap.f[, c("Mcap.Ct.sd", "C.Ct.sd", "D.Ct.sd")], horizontal=T, las=2, ylim=c(0,5),
#         main="Filtered runs")

# Identify overall dominant symbiont clade across time points
syms <- ldply(split(Mcap.f[, "dom"], f = Mcap.f$sample), table)
dom <- apply(syms, 1, FUN=function(x) ifelse(x[2] > x[3], "C", ifelse(x[3] > x[2], "D", "CD")))
dom <- cbind(syms, dom)
rownames(dom) <- dom$".id"
Mcap.f$tdom <- dom[as.character(Mcap.f$sample), "dom"]
# -------------------------------------------------------------------------------------------------
# • Filter low quality data for quantitative analyses (samples with delayed host amplification) -----------------
# # Visualize distribution of Mcap Ct values
# par(mfrow=c(2,1), mar=c(0,3,5,1))
# hist(Mcap.f$Mcap.Ct.mean, breaks=seq(15,40,1), xlim=c(15,40), xaxt="n", main="Mcap Ct values")
# par(mar=c(5,3,0,1))
# boxplot(Mcap.f$Mcap.Ct.mean,ylab="Ct value", horizontal=T, ylim=c(15,40), frame=F)
# Remove samples with high outlying Mcap Ct values for quantitative analyses
thresh <- boxplot.stats(Mcap.f$Mcap.Ct.mean)$stats[5]
Mcap.ff <- Mcap.f[which(Mcap.f$Mcap.Ct.mean <= thresh), ]
# -------------------------------------------------------------------------------------------------
