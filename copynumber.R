# Import data
library(reshape2)
source("~/Documents/Academia/HIMB/steponeR/steponeR.R")
files <- list("20150911_KBayRecov_Mcap_copyno_1_data.csv", 
              "20150911_KBayRecov_Mcap_copyno_2_data_mod.csv")
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

# Import cell count data
cells <- read.csv("cells.csv")
cells$cells.µL <- cells$extract.cells * 0.955 / 50  # Calculate cells per µL given 95.5% extraction efficiency and 50 µL elution volume

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
copynumber




