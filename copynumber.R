# Import data
file <- "20150907_KBayRecov_Mcap_copynumber_data.csv"
cn <- read.csv(file, skip=6, na.strings="Undetermined")[, c(2:4, 7, 10)]
cn <- split(cn, f=cn$Task)
cn <- lapply(cn, droplevels)

# Plot standard curves with Quantity as response variable
plot(log10(Quantity) ~ C_, data=cn$STANDARD, col=c("blue", "red")[cn$STANDARD$Target.Name])

# Fit standard curve model
std <- lm(log10(Quantity) ~ C_ + Target.Name, data=cn$STANDARD)
anova(std)  # No difference between Targets (b/c different automatic CT thresholds for C and D)
std <- lm(log10(Quantity) ~ C_, data=cn$STANDARD)
# Remove outliers with residuals > 2 s.d.'s from fit
out <- romr.fnc(std, data=cn$STANDARD, trim=2)
std <- update(std, data=out$data)
summary(std)

# Plot data without outliers
stddat <- model.frame(std)
plot(`log10(Quantity)` ~ C_, data=stddat, col=c("blue", "red")[cn$STANDARD$Target.Name])

# Plot fitted values
stdcurve <- data.frame(C_=seq(15,35,1))
stdcurve$fit <- predict(std, stdcurve)
with(stdcurve, lines(C_, fit))

# Calculate unknown quantities based on standard curve
cn$UNKNOWN$fit <- 10^predict(std, cn$UNKNOWN)
cn$UNKNOWN
hist(cn$UNKNOWN$C_)

# Restructure data
cndat <- dcast(cn$UNKNOWN, formula = Sample.Name ~ Target.Name, value.var="Quantity", fun.aggregate = mean)
cndat$Sample <- substr(cndat$Sample.Name, 5, 7)
cndat

# Import cell count data
cells <- read.csv("cells.csv")
cells$cells <- cells$cells.mL * 0.5 * (100/1000) * 0.955 * (1/50)  # Actual # of cells going into reaction

final <- merge(cndat, cells, by="Sample")
final[is.na(final)] <- 0

#D <- final[which(final$C==0), ]
#final <- D
#mean(final$D / final$cells)

# Singular value decomposition
copies<-as.matrix(final[,3:4])
copies[is.na(copies)] <- 0
copies
cells<-as.matrix(final[,6])
cells

svd<-svd(copies)
svd
wp<-matrix(c(1/svd$d[1],0,0,1/svd$d[2]),ncol=2)
wp

x=svd$v%*%wp%*%t(svd$u)%*%cells

copynumber=1/x
copynumber

#----------
# Restructure data
cndat <- dcast(cn$UNKNOWN, formula = Sample.Name ~ Target.Name, value.var="C_", fun.aggregate = mean)
cndat$Sample <- substr(cndat$Sample.Name, 5, 7)
cndat

cells <- read.csv("cells.csv")
cells$cells <- cells$cells.mL * 0.5 * (100/1000) * 0.955 * (1/50)  # Actual # of cells going into reaction

final <- merge(cndat, cells, by="Sample")
final[is.na(final)] <- 0
final

D <- final[which(final$D < final$C), ]
final <- D
final <- final[-c(2, 4), ]

plot(D ~ log(cells, 2), final)
