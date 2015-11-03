# Title: Patterns of bleaching and recovery of Montipora capitata in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 26 October, 2015
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
