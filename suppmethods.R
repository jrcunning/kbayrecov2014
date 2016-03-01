# Title: Patterns of bleaching and recovery of Montipora capitata in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 1 March 2016
# =================================================================================================
# Supplementary Materials
# =================================================================================================

# Fluorescence normalization --------------------
file <- "data/supp/20150907_KBayRecov_Mcap_fluornorm_data.csv"
fn <- read.csv(file, skip=6, na.strings="Undetermined")[, c(1:3, 7, 10)]
fn <- droplevels(fn[which(fn$Sample.Name!=""), ])

# Plot standard curves
plot(C_ ~ log10(Quantity), data=fn, pch=21, bg=fn$Target.Name, ylab="CT value")
legend("topright", levels(fn$Target.Name), pch=21, pt.bg=factor(levels(fn$Target.Name)))

# Set clade D target as baseline contrast
contrasts(fn$Target.Name) <- contr.treatment(levels(fn$Target.Name), base=2)
# Fit model with predictors quantity and target
mod <- lm(C_ ~ log10(Quantity) * Target.Name, data=fn) 
anova(mod)  # No interaction = slopes equal
mod <- update(mod, ~ log10(Quantity) + Target.Name) # Remove interaction
anova(mod)
summary(mod)
coef(mod)  # Fluorescence normalization values: 0 for clade D, 0.84815 for Mcap, and 2.26827 for clade C

# Plot fitted values for each target
newdat <- expand.grid(Quantity=10^seq(2,6,1), Target.Name=factor(c("D", "C", "Mcap")))
newdat$fit <- predict(mod, newdat)
with(newdat[which(newdat$Target.Name=="D"), ], lines(log10(Quantity), fit, col="red"))
with(newdat[which(newdat$Target.Name=="Mcap"), ], lines(log10(Quantity), fit, col="green"))
with(newdat[which(newdat$Target.Name=="C"), ], lines(log10(Quantity), fit, col="black"))

# Copy number -----------------
# Import data
library(reshape2)
source("~/Documents/Academia/HIMB/steponeR/steponeR.R")
files <- list("data/supp/20150911_KBayRecov_Mcap_copyno_1_data.csv", 
              "data/supp/20150911_KBayRecov_Mcap_copyno_2_data.csv")
cn <- steponeR(files)  # Imports data, calculates copies in unknowns based on standards for each plate

# Plot standards
with(cn$standards$data, {
  plot(CT ~ log10(Quantity),
       col=Target.Name, pch=seq(1:nlevels(Filename))[Filename])
  legend("topright", legend=levels(interaction(Target.Name, Filename)), 
         col=levels(Target.Name), pch=seq(1:nlevels(Filename))[levels(Filename)], cex=0.5)
})

# Cast data frame averaging copies among technical replicates
cndat <- dcast(cn$unknowns, Sample.Name + Filename ~ Target.Name, mean, na.rm=T, value.var="copies")
# Separate sample names into colony and extraction replicate
cndat <- data.frame(colsplit(cndat$Sample.Name, "-", names=c("colony", "extract.rep")), cndat)

# Plate 1 had 1 µL template loaded, plate 2 had 2 µL template loaded
cndat$template.µL <- ifelse(cndat$Filename=="20150911_KBayRecov_Mcap_copyno_1_data.csv", 1, 2)

# Cell count data: number of cells from which DNA was extracted for each sample
cells <- data.frame(colony=c(52,72,58,80,57,69),
                    extract.cells=c(100000,100000,83295,71595,94545,100000))
cells$cells.µL <- cells$extract.cells * 0.955 / 50  # Calculate cells per µL in extracted DNA given 95.5% extraction efficiency and 50 µL elution volume

# Merge copies and cells data, calculate expected cells per reaction
cndat <- merge(cndat, cells[, c("colony", "cells.µL")], by="colony")
cndat$cells.rxn <- cndat$cells.µL * cndat$template.µL
cndat[is.na(cndat)] <- 0


# Singular value decomposition and least squares solution
copies <- as.matrix(cndat[, c("C", "D")])
cells <- as.matrix(cndat[,"cells.rxn"])
svd <- svd(copies)
wp <- matrix(c(1 / svd$d[1], 0, 0, 1/svd$d[2]), ncol=2)
x <- svd$v %*% wp %*% t(svd$u) %*% cells
copynumber <- 1 / x
copynumber  # Values of 33 and 3 utilized for data analysis

# ITS2 sequence analysis ----------------

library(phyloseq)

# Import data as phyloseq object
OTU97 <- otu_table(read.table("data/ITS2/OTUs_97/otu_table_97.tsv", header=T, check.names=F, 
                              row.names=1, skip=1, comment.char="", sep="\t"), taxa_are_rows=T)
OTU100 <- otu_table(read.table("data/ITS2/OTUs_100/otu_table_100.tsv", header=T, check.names=F, 
                               row.names=1, skip=1, comment.char="", sep="\t"), taxa_are_rows=T)
TAX97 <- read.table("data/ITS2/OTUs_97/blast_taxonomy_97/rep_set_97_tax_assignments.txt", 
                    sep="\t", row.names=1, col.names=c("OTU.ID", "Taxonomy", "Eval", "Subtype"))
TAX97$otuname <- rownames(TAX97)
TAX97 <- tax_table(as.matrix(TAX97))
TAX100 <- read.table("data/ITS2/OTUs_100/blast_taxonomy_100/rep_set_100_tax_assignments.txt", sep="\t", row.names=1,
                     col.names=c("OTU.ID", "Taxonomy", "Eval", "Subtype"))
TAX100$otuname <- rownames(TAX100)
TAX100 <- tax_table(as.matrix(TAX100))
SAM <- sample_data(read.csv("data/ITS2/sample_data.csv", row.names=1))
phy97 <- phyloseq(OTU97, TAX97, SAM)
phy100 <- phyloseq(OTU100, TAX100, SAM)
phy97; phy100

# Remove taxa not seen at least 3 times in at least one sample
phy97.f <- filter_taxa(phy97, function(x) sum(x >= 3) > 1, prune=T)
phy100.f <- filter_taxa(phy100, function(x) sum(x >= 3) > 1, prune=T)
phy97.f; phy100.f
# Transform sample counts to relative abundance
phy97.f.p <- transform_sample_counts(phy97.f, function(x) x/sum(x))
phy100.f.p <- transform_sample_counts(phy100.f, function(x) x/sum(x))
# Barplots
plot_bar(phy97.f.p, fill="Subtype", facet_grid=~vis)
plot_bar(phy100.f.p, fill="otuname", facet_grid=~vis)

# Most abundant OTUs
head(otu_table(phy97.f.p)[order(rowSums(otu_table(phy97.f.p)), decreasing=T)])
tax_table(phy97.f.p)[rownames(head(otu_table(phy97.f.p)[order(rowSums(otu_table(phy97.f.p)), decreasing=T)]))]
head(otu_table(phy100.f.p)[order(rowSums(otu_table(phy100.f.p)), decreasing=T)])
tax_table(phy100.f.p)[rownames(head(otu_table(phy100.f.p)[order(rowSums(otu_table(phy100.f.p)), decreasing=T)]))]

# Comparisons using PERMANOVA
library(vegan)
# Compare C-dominated vs. D-dominated colonies: 97% OTUs
perm.97.clade <- adonis(phyloseq::distance(phy97.f.p, "bray") ~ sym, 
                        data=as(sample_data(phy97.f.p), "data.frame"), permutations=999)
# Compare bleached C-dominated colonies vs. non-bleached C-dominated colonies: 97% OTUs
phy97.f.p.C <- subset_samples(phy97.f.p, sym=="C")
perm.97.visC <- adonis(phyloseq::distance(phy97.f.p.C, "bray") ~ vis, 
                       data=as(sample_data(phy97.f.p.C), "data.frame"), permutations=999)
# Compare C-dominated vs. D-dominated colonies: 100% OTUs
perm.100.clade <- adonis(phyloseq::distance(phy100.f.p, "bray") ~ sym, 
                        data=as(sample_data(phy100.f.p), "data.frame"), permutations=999)
# Compare bleached C-dominated colonies vs. non-bleached C-dominated colonies: 100% OTUs
phy100.f.p.C <- subset_samples(phy100.f.p, sym=="C")
perm.100.visC <- adonis(phyloseq::distance(phy100.f.p.C, "bray") ~ vis, 
                       data=as(sample_data(phy100.f.p.C), "data.frame"), permutations=999)
# Compare non-bleached C-dominated colonies vs. non-bleached D-dominated colonies: 97% OTUs
phy97.f.p.NB <- subset_samples(phy97.f.p, vis=="not bleached")
perm.97.cladeNB <- adonis(phyloseq::distance(phy97.f.p.NB, "bray") ~ sym, 
                      data=as(sample_data(phy97.f.p.NB), "data.frame"), permutations=999)
# Compare non-bleached C-dominated colonies vs. non-bleached D-dominated colonies: 100% OTUs
phy100.f.p.NB <- subset_samples(phy100.f.p, vis=="not bleached")
perm.100.cladeNB <- adonis(phyloseq::distance(phy100.f.p.NB, "bray") ~ sym, 
                          data=as(sample_data(phy100.f.p.NB), "data.frame"), permutations=999)

# Figure S1
pdf(file = "output/FigureS1.pdf", width = 6.65354, height=7.65354)
par(mfrow=c(2,1), mar=c(2,3,4,1), mgp=c(3,0.5,0))

# 97%-OTUs barplot
samdat <- data.frame(sample_data(phy97.f.p))
samdat <- samdat[c("6K","20K","58K","66K","110K","128K","3K","11K","127K"), ]
typerelabund <- as.matrix(otu_table(phy97.f.p)[order(data.frame(tax_table(phy97.f.p))$Subtype), rownames(samdat)])
# Get info for plotting OTU names and blast hits on top of barplot
blasthits <- as.character(data.frame(tax_table(phy97.f.p))[order(data.frame(tax_table(phy97.f.p))$Subtype), "Subtype"])
names <- as.character(data.frame(tax_table(phy97.f.p))[order(data.frame(tax_table(phy97.f.p))$Subtype), "otuname"])
divides <- apply(typerelabund, 2, cumsum)
heights <- diff(divides)
heights[heights < 0.04] <- NA
bars <- barplot(typerelabund, col=gray.colors(10), las=1, xlab="", ylab="", cex.axis=0.75, cex.names=0.75)
mtext(side=2, text = "Relative Abundance", line=2)
for (i in 1:length(bars)) {
  text(rep(bars[i], length(heights[which(!is.na(heights[,i])),i])), 
       divides[which(!is.na(heights[,i])),i] + heights[which(!is.na(heights[,i])),i] / 2, 
       labels=paste(names[which(!is.na(heights[,i]))+1], blasthits[which(!is.na(heights[,i]))+1],sep="\n"),
       cex=0.4)
}
text(bars, rep(par("usr")[4], length(bars)),
     labels=c("D","D","D","C","C","C","C","C","C"), xpd=T, pos=3, cex=0.75)
text(bars, rep(par("usr")[4]+0.075, length(bars)),
     labels=c("NB","NB","NB","NB","NB","NB","B","B","B"), xpd=T, pos=3, cex=0.75)
mtext(side=3, "A.) 97% OTUs", xpd=T, adj=0, line=2.5, font=2)

# 100%-OTUs barplot 
par(mar=c(2,3,4,1))
samdat <- data.frame(sample_data(phy100.f.p))
samdat <- samdat[c("6K","20K","58K","66K","110K","128K","3K","11K","127K"), ]
typerelabund <- as.matrix(otu_table(phy100.f.p)[order(data.frame(tax_table(phy100.f.p))$Subtype), rownames(samdat)])

# Get info for plotting OTU names and blast hits on top of barplot
blasthits <- as.character(data.frame(tax_table(phy100.f.p))[order(data.frame(tax_table(phy100.f.p))$Subtype), "Subtype"])
names <- as.character(data.frame(tax_table(phy100.f.p))[order(data.frame(tax_table(phy100.f.p))$Subtype), "otuname"])
divides <- apply(typerelabund, 2, cumsum)
heights <- diff(divides)
heights[heights < 0.04] <- NA

# Plot barplot and OTU names and blast hits
bars <- barplot(typerelabund, col=rainbow(941), border=NA, las=1, xlab="Sample", ylab="", cex.axis=0.75, cex.names=0.75)
mtext(side=2, text = "Relative Abundance", line=2)
for (i in 1:length(bars)) {
  text(rep(bars[i], length(heights[which(!is.na(heights[,i])),i])), 
       divides[which(!is.na(heights[,i])),i] + heights[which(!is.na(heights[,i])),i] / 2, 
       labels=paste(names[which(!is.na(heights[,i]))+1],blasthits[which(!is.na(heights[,i]))+1],sep="\n"),
       cex=0.4)
}
text(bars, rep(par("usr")[4], length(bars)),
     labels=c("D","D","D","C","C","C","C","C","C"), xpd=T, pos=3, cex=0.75)
text(bars, rep(par("usr")[4]+0.075, length(bars)),
     labels=c("NB","NB","NB","NB","NB","NB","B","B","B"), xpd=T, pos=3, cex=0.75)
mtext(side=3, "B.) 100% OTUs", xpd=T, adj=0, line=2.5, font=2)
dev.off()

# Temperature and light data analysis ------
library(scales); library(zoo)
reefcols <- c("#8dd3c7", "#bebada", "#d9d9d9")
 
# Import temperature data
rf25.temp <- rbind(read.csv("data/temp_light/temp_reef25_2014.csv"),
                   read.csv("data/temp_light/temp_reef25_2015.csv"))
rf44.temp <- rbind(read.csv("data/temp_light/temp_reef44_2014.csv"),
                   read.csv("data/temp_light/temp_reef44_2015.csv"))
HIMB.temp <- rbind(read.csv("data/temp_light/temp_reefHIMB_2014.csv"),
                   read.csv("data/temp_light/temp_reefHIMB_2015.csv"))
# set date to date class
rf25.temp$date <- as.Date(rf25.temp$date, format="%m/%e/%Y")
rf44.temp$date <- as.Date(rf44.temp$date, format="%m/%e/%Y")
HIMB.temp$date <- as.Date(HIMB.temp$date, format="%m/%e/%Y")
# Subset data only up until threshdate
threshdate <- as.Date("2015-05-06", format="%F")
rf25.temp <- rf25.temp[which(rf25.temp$date < threshdate), ]
rf44.temp <- rf44.temp[which(rf44.temp$date < threshdate), ]
HIMB.temp <- HIMB.temp[which(HIMB.temp$date < threshdate), ]

# Aggregate temperature data by daily mean, minimum, and maximum
rf25t.split <- split(rf25.temp, f=rf25.temp$date < as.Date("2014-11-21", format="%F"))
rf25t.1 <- aggregate(data.frame(mean=rf25t.split[[2]]$temp_raw_C), by=list(date=rf25t.split[[2]]$date), FUN=mean)
rf25t.1$min <- aggregate(rf25t.split[[2]]$temp_raw_C, by=list(rf25t.split[[2]]$date), FUN=min)$x
rf25t.1$max <- aggregate(rf25t.split[[2]]$temp_raw_C, by=list(rf25t.split[[2]]$date), FUN=max)$x
rf25t.2 <- aggregate(data.frame(mean=rf25t.split[[1]]$temp_cal_C), by=list(date=rf25t.split[[1]]$date), FUN=mean)
rf25t.2$min <- aggregate(rf25t.split[[1]]$temp_cal_C, by=list(rf25t.split[[1]]$date), FUN=min)$x
rf25t.2$max <- aggregate(rf25t.split[[1]]$temp_cal_C, by=list(rf25t.split[[1]]$date), FUN=max)$x
rf44t <- aggregate(data.frame(mean=rf44.temp$temp_cal), by=list(date=rf44.temp$date), FUN=mean)
rf44t$min <- aggregate(rf44.temp$temp_cal, by=list(rf44.temp$date), FUN=min)$x
rf44t$max <- aggregate(rf44.temp$temp_cal, by=list(rf44.temp$date), FUN=max)$x
HIMBt <- aggregate(data.frame(mean=HIMB.temp$temp_cal), by=list(date=HIMB.temp$date), FUN=mean)
HIMBt$min <- aggregate(HIMB.temp$temp_cal, by=list(HIMB.temp$date), FUN=min)$x
HIMBt$max <- aggregate(HIMB.temp$temp_cal, by=list(HIMB.temp$date), FUN=max)$x

# Import light data
rf25.light <- read.csv("data/temp_light/light_reef25_2015.csv")
rf44.light <- rbind(read.csv("data/temp_light/light_reef44_2014.csv"),
                    read.csv("data/temp_light/light_reef44_2015.csv"))
HIMB.light <- rbind(read.csv("data/temp_light/light_reefHIMB_2014.csv"),
                    read.csv("data/temp_light/light_reefHIMB_2015.csv"))
# Set date to date class
rf25.light$date <- as.Date(rf25.light$date, format="%e/%m/%Y")
rf44.light$date <- as.Date(rf44.light$date, format="%e/%m/%Y")
HIMB.light$date <- as.Date(HIMB.light$date, format="%e/%m/%Y")
# Subset data only up until threshdate
threshdate <- as.Date("2015-08-06", format="%F")
rf25.light <- rf25.light[which(rf25.light$date < threshdate), ]
rf44.light <- rf44.light[which(rf44.light$date < threshdate), ]
HIMB.light <- HIMB.light[which(HIMB.light$date < threshdate), ]
# Aggregate light data by daily mean, minimum, and maximum
rf25l <- aggregate(data.frame(mean=rf25.light$PAR_calibrated_umol.m.2.s), by=list(date=rf25.light$date), FUN=mean)
rf25l$max <- aggregate(rf25.light$PAR_calibrated_umol.m.2.s, by=list(date=rf25.light$date), FUN=max)$x
rf44l <- aggregate(data.frame(mean=rf44.light$PAR_calibrated_umol.m.2.s), by=list(date=rf44.light$date), FUN=mean)
rf44l$max <- aggregate(rf44.light$PAR_calibrated_umol.m.2.s, by=list(date=rf44.light$date), FUN=max)$x
HIMBl <- aggregate(data.frame(mean=HIMB.light$PAR_calibrated_umol.m.2.s), by=list(date=HIMB.light$date), FUN=mean)
HIMBl$max <- aggregate(HIMB.light$PAR_calibrated_umol.m.2.s, by=list(date=HIMB.light$date), FUN=max)$x
rf25l$dli <- rf25l$mean * 0.0864  # Convert from µmol/s to mol/d
rf44l$dli <- rf44l$mean * 0.0864
HIMBl$dli <- HIMBl$mean * 0.0864

# # Plot daily light integral for each reef
# k=5 # k-day moving averages
# plot(dli ~ date, rf44l, type="n", main="Temperature mean and range", ylab="°C")
# legend("topright", lty=1, col=c(reefcols[1:2], "darkgray"), legend=c("44","25","HIMB"), lwd=2)
# with(na.omit(data.frame(date=HIMBl$date, dli=rollmean(HIMBl$dli, k, fill=NA))), lines(date, dli, col=reefcols[3], lwd=1))
# with(na.omit(data.frame(date=rf25l$date, dli=rollmean(rf25l$dli, k, fill=NA))), lines(date, dli, col=reefcols[2], lwd=1))
# with(na.omit(data.frame(date=rf44l$date, dli=rollmean(rf44l$dli, k, fill=NA))), lines(date, dli, col=reefcols[1], lwd=1))


# Recalibrate reef 25 using linear fit
# equation from excel: y=0.1178x
rf25.light$PAR_linear_calibrated_umol.m.2.s <- rf25.light$PAR_raw * 0.1178
rf25l.lin <- aggregate(data.frame(mean=rf25.light$PAR_linear_calibrated_umol.m.2.s), by=list(date=rf25.light$date), FUN=mean)
rf25l.lin$max <- aggregate(rf25.light$PAR_linear_calibrated_umol.m.2.s, by=list(date=rf25.light$date), FUN=max)$x
rf25l.lin$dli <- rf25l.lin$mean * 0.0864 

# Recalibrate reef 44 using linear fit
# equation from excel: y=0.0602x
rf44.light$PAR_linear_calibrated_umol.m.2.s <- rf44.light$PAR_raw * 0.0602
rf44l.lin <- aggregate(data.frame(mean=rf44.light$PAR_linear_calibrated_umol.m.2.s), by=list(date=rf44.light$date), FUN=mean)
rf44l.lin$max <- aggregate(rf44.light$PAR_linear_calibrated_umol.m.2.s, by=list(date=rf44.light$date), FUN=max)$x
rf44l.lin$dli <- rf44.lin$mean * 0.0864 

# Recalibrate reef HIMB using linear fit
# logger 2489 (17 Dec 2014 - 26 May 2015): << calibration pending >>
# logger 4805 (26 May 2015 - 31 Dec 2015): y=0.0709x
HIMB.light[which(HIMB.light$date < as.Date("2015-05-26")), "PAR_linear_calibrated_umol.m.2.s"] <- NA #HIMB.light[which(HIMB.light$date < as.Date("2015-05-26")), "PAR_raw"] * ______
HIMB.light[which(HIMB.light$date >= as.Date("2015-05-26")), "PAR_linear_calibrated_umol.m.2.s"] <- HIMB.light[which(HIMB.light$date >= as.Date("2015-05-26")), "PAR_raw"] * 0.0709
HIMBl.lin <- aggregate(data.frame(mean=HIMB.light$PAR_linear_calibrated_umol.m.2.s), by=list(date=HIMB.light$date), FUN=mean, na.rm=T)
HIMBl.lin$max <- aggregate(HIMB.light$PAR_linear_calibrated_umol.m.2.s, by=list(date=HIMB.light$date), FUN=max)$x
HIMBl.lin$dli <- HIMBl.lin$mean * 0.0864 

# Plot temp. daily mean, min, and max for each reef
pdf(file="output/FigureS2.pdf",  width = 6.65354, height=5)
par(mfrow=c(2,1), mar=c(2,3,1,1), mgp=c(2,0.5,0))
k=1; lwd=0.5 # k-day moving averages
plot(mean ~ date, rf44t, type="n", ylab="Temperature (°C)", xaxt="n", xlab="")
mtext(expression(bold("A")), 2, adj=4.5, las=1, padj=-8)
axis.Date(1, at=seq(min(rf44t$date), max(rf44t$date), by="1 mon"), format="%b '%y")
legend("topright", lty=1, col=c(reefcols[1:2], "darkgray"), legend=c("44","25","HIMB"), lwd=2, bty="n")
with(na.omit(data.frame(date=HIMBt$date, rollmean(HIMBt[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[3], 0.3))
  lines(date, mean, col=reefcols[3], lwd=lwd)
})
with(na.omit(data.frame(date=rf25t.1$date, rollmean(rf25t.1[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[2], 0.3))
  lines(date, mean, col=reefcols[2], lwd=lwd)
})
with(na.omit(data.frame(date=rf25t.2$date, rollmean(rf25t.2[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[2], 0.3))
  lines(date, mean, col=reefcols[2], lwd=lwd)
})
with(na.omit(data.frame(date=rf44t$date, rollmean(rf44t[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[1], 0.3))
  lines(date, mean, col=reefcols[1], lwd=lwd)
})
k=7; lwd=2
with(na.omit(data.frame(date=HIMBt$date, rollmean(HIMBt[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[3], 0.3))
  lines(date, mean, col=reefcols[3], lwd=lwd)
})
with(na.omit(data.frame(date=rf25t.1$date, rollmean(rf25t.1[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[2], 0.3))
  lines(date, mean, col=reefcols[2], lwd=lwd)
})
with(na.omit(data.frame(date=rf25t.2$date, rollmean(rf25t.2[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[2], 0.3))
  lines(date, mean, col=reefcols[2], lwd=lwd)
})
with(na.omit(data.frame(date=rf44t$date, rollmean(rf44t[, c("mean", "min", "max")], k, fill=NA))), {
  #polygon(x=c(date, rev(date)), y=c(max, rev(min)), border=NA, col=alpha(reefcols[1], 0.3))
  lines(date, mean, col=reefcols[1], lwd=lwd)
})
# Plot LINEARLY CALIBRATED daily light integral for each reef
k=1; lwd=0.5 # k-day moving averages
plot(dli ~ date, rf44l.lin, type="n", ylab="DLI (mol m-2 d-1)", xaxt="n", xlab="")
mtext(expression(bold("B")), 2, adj=4.5, las=1, padj=-8)
axis.Date(1, at=seq(min(rf44l.lin$date), max(rf44l.lin$date), by="1 mon"), format="%b '%y")
#mtext(side=3, text="HIMB - logger 4805 only")
legend("topleft", lty=1, col=c(reefcols[1:2], "darkgray"), legend=c("44","25","HIMB"), lwd=2, bty="n")
with(na.omit(data.frame(date=HIMBl.lin[which(HIMBl.lin$date >= as.Date("2015-05-26")), "date"], 
                        dli=rollmean(HIMBl.lin[which(HIMBl.lin$date >= as.Date("2015-05-26")), "dli"], k, fill=NA))), 
     lines(date, dli, col=reefcols[3], lwd=lwd))
with(na.omit(data.frame(date=rf25l.lin$date, dli=rollmean(rf25l.lin$dli, k, fill=NA))), lines(date, dli, col=reefcols[2], lwd=lwd))
with(na.omit(data.frame(date=rf44l.lin$date, dli=rollmean(rf44l.lin$dli, k, fill=NA))), lines(date, dli, col=reefcols[1], lwd=lwd))
k=7; lwd=2
with(na.omit(data.frame(date=HIMBl.lin[which(HIMBl.lin$date >= as.Date("2015-05-26")), "date"], 
                        dli=rollmean(HIMBl.lin[which(HIMBl.lin$date >= as.Date("2015-05-26")), "dli"], k, fill=NA))), 
     lines(date, dli, col=reefcols[3], lwd=lwd))
with(na.omit(data.frame(date=rf25l.lin$date, dli=rollmean(rf25l.lin$dli, k, fill=NA))), lines(date, dli, col=reefcols[2], lwd=lwd))
with(na.omit(data.frame(date=rf44l.lin$date, dli=rollmean(rf44l.lin$dli, k, fill=NA))), lines(date, dli, col=reefcols[1], lwd=lwd))
dev.off()
