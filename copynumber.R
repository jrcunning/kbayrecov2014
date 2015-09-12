# Import data
#file <- "20150907_KBayRecov_Mcap_copynumber_data.csv"
file <- "20150911_KBayRecov_Mcap_copyno_1_data.csv"
cn <- read.csv(file, skip=6, na.strings="Undetermined")[, c(2:4, 7, 10)]
cn <- split(cn, f=cn$Task)
cn <- lapply(cn, droplevels)

# Restructure data #################
library(reshape2)
cndat <- dcast(cn$UNKNOWN, formula = Sample.Name ~ Target.Name, value.var="Quantity", fun.aggregate = mean)
cndat <- data.frame(colsplit(cndat$Sample.Name, "-", names=c("colony", "extract.rep")), cndat[, -1])
cndat

# Import cell count data
cells <- read.csv("cells.csv")
cells$cells.µL <- cells$extract.cells * 0.955 /50  # Calculate cells per µL given 95.5% extraction efficiency and 50 µL elution volume
cells$cells.µL

# Merge copies and cells data
cndat <- merge(cndat, cells, by="colony")
cndat[is.na(cndat)] <- 0



# Singular value decomposition
copies<-as.matrix(cndat[,3:4])
copies[is.na(copies)] <- 0
copies
cells<-as.matrix(cndat[,6])
cells

svd<-svd(copies)
svd
wp<-matrix(c(1/svd$d[1],0,0,1/svd$d[2]),ncol=2)
wp

x=svd$v%*%wp%*%t(svd$u)%*%cells

copynumber=1/x
copynumber


