# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 27 August, 2015

# Run setup script
source("setup.R")
# Set seed
set.seed(39059978)

# =================================================================================================
# SYMBIONT COMMUNITY STRUCTURE AND BLEACHING
# =================================================================================================
# • Analysis: Overall symbiont clade composition in samples and colonies --------------------------
# Presence of clade C and D in individual samples
symtab <- table(Mcap.f$syms)
samples <- c(symtab[1], symtab[2] + symtab[3], symtab[4])
# Presence of clade C and D in colonies (detected at least once)
colonies <- aggregate(Mcap.f$syms, by=list(colony=Mcap.f$colony), FUN=paste, collapse="")
colonies$C[grep("C", colonies$x)] <- "C"
colonies$D[grep("D", colonies$x)] <- "D"
colonies$present <- ifelse(is.na(colonies$C), ifelse(is.na(colonies$D), "none", "D only"), ifelse(is.na(colonies$D), "C only", "C+D"))
# Summary of clade composition of samples and colonies
clades <- data.frame(Colonies=matrix(prop.table(table(colonies$present))), Samples=prop.table(samples))
# Proportion of colonies with C or D dominant over time
propdom <- prop.table(table(Mcap.f[which(Mcap.f$fdate=="20141024"), "tdom"]))
# Proportion D in mixed samples
propD <- Mcap.f$propD[which(Mcap.f$propD > 0 & Mcap.f$propD < 1)]
range(propD)
# Percent of samples with >10% non-dominant symbiont (between 10% and 90% clade D)
sum(prop.table(hist(propD, plot=F)$counts)[2:9])

# These tests are not relevant because population was not randomly sampled (targeted bleached and healthy pairs)
# # Analyze patterns in clade composition across reefs
# cladecomp <- merge(unique(Mcap.f[, c("colony", "reef", "tdom")]), colonies)
# # Analyze proportions of mixed colonies across reefs
# chisq1 <- chisq.test(cladecomp$reef, cladecomp$present)  # no differences
# chisq1$observed; chisq1$p.value  # No differences
# # Analyze proportions of dominant symbionts across reefs
# chisq2 <- chisq.test(cladecomp$reef, cladecomp$tdom)
# chisq2$observed; chisq2$p.value  # No differences
# • Figure 1: Overall symbiont community composition in Montipora capitata  ---------------
pdf(file="Figure1.pdf", height=3, width=3)
# Plot histogram of proportion clade D in mixed C+D samples
par(mfrow=c(1,2), mar=c(10,1.5,1,2), mgp=c(1,0.25,0), tcl=-0.25)
hist(propD, plot=T, main="", xaxt="n", xlab="", yaxt="n", ylab="", col = "gray60")
box()
axis(side=1, at=seq(0,1,0.1), labels=NA, tck=-0.05)
axis(side=1, at=c(0,0.5,1), labels=c(">0","0.5","<1"), cex.axis=0.6, mgp=c(1,0,0), tck=0)
mtext(side=1, "Prop. clade D", line=1, cex=0.75)
axis(side=2, tck=-0.05, cex.axis=0.6, mgp=c(1,0.15,0))
mtext(side=2, "# samples", line=0.8, cex=0.75)
# Save coordinates from histogram plot for drawing lines
save1.x <- grconvertX(par("usr")[2], from='user', to='ndc')
save1.y <- grconvertY(par("usr")[3], from='user', to='ndc')
save2.y <- grconvertY(par("usr")[4], from='user', to='ndc')
# Add legend for barplot
legend("bottom", inset=c(0,-2), legend=c("D only", "C and D", "C only"), xpd=NA, cex=1,
       fill=c("gray95", "gray60", "gray20"), bty="n")
# Plot barplot of clade composition of samples and colonies
par(mar=c(2,0,1,2), mgp=c(1.5,0.25,0), lwd=1)
b1 <- barplot(as.matrix(clades)[,c(2,1)], beside=F, space=0.75, mgp=c(1,0.25,0),
              col = c("gray20", "gray60", "gray95"), axes=F, cex.names=0.75, line=-0.3, names.arg=c("",""))
axis(side=1, at=b1, labels=c("Samples", "Colonies"), lty=0, cex.axis=0.75, padj=1, line=-0.6)
text(b1, par("usr")[4] - 0.03, pos=3, xpd=T, cex=0.75,
     labels = paste("n=", c(sum(samples), dim(colonies)[1]), sep=""))
text(b1[2], par("usr")[3] - 0.11, labels = "(sampled 3-6x)", cex=0.5, xpd=T)
axis(side=2, pos=quantile(par("usr")[1:2], 0.54), las=1, mgp=c(0,-0.4,0), lty=0, cex.axis=0.6,
     at=c(0,0.2,0.4,0.6,0.8,1),
     labels = c("-0.0-", "-0.2-", "-0.4-", "-0.6-", "-0.8-", "-1.0-"))
segments(x0 = par("usr")[1], y0 = clades[1,2], 
         x1 = grconvertX(save1.x, from='ndc'), 
         y1 = grconvertY(save1.y, from='ndc'), lty=1, xpd=NA)
segments(x0 = par("usr")[1], y0 = 1 - clades[3,2], 
         x1 = grconvertX(save1.x, from='ndc'), 
         y1 = grconvertY(save2.y, from='ndc'), lty=1, xpd=NA)
brackets(x1=par("usr")[2], y1=propdom["C"], x2=par("usr")[2], y2=par("usr")[3],
         h=0.3, xpd=NA, type=1)
brackets(x1=par("usr")[2], y1=par("usr")[4], x2=par("usr")[2], y2=propdom["C"],
         h=0.3, xpd=NA, type=1)
text(x=rep(par("usr")[2] + 0.5, 2), y=c(0.4, 0.9), labels = c("C dominant", "D dominant"), 
     xpd=T, pos=1, srt=90, cex=0.75)
dev.off()
# • Analysis: Relationship between dominant symbiont clade and bleaching -------------------------------
# Summarize data - reef, visual status, and dominant symbiont for each colony
Mcap.f.summ <- unique(Mcap.f[, c("colony", "vis", "reef", "tdom")])
# Test for differences in dominant clade between bleached and unbleached colonies
vis.chi <- chisq.test(Mcap.f.summ$vis, Mcap.f.summ$tdom)
vis.chi[c("observed", "p.value")]

# Test for differences in dominant clade of unbleached colonies across reefs
with(Mcap.f.summ[which(Mcap.f.summ$vis=="not bleached"), ], {
  chisq.test(reef, tdom)[c("observed", "p.value")]
})  # No difference

# Analyze bleaching severity (S/H cell ratio) in October (peak of bleaching)
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
bleach <- lm(log(tot.SH) ~ tdom:vis, data=Mcap.ff.oct)
anova(bleach, test="F")  # tdom:vis significant, (reef is not significant if included)
TukeyHSD(aov(bleach))  # Pairwise tests between groups
eff <- data.frame(effect("tdom:vis", bleach), confidence.level=0.95)[c(1,3,4),]
# Means and standard errors by vis only
means <- aggregate(log(Mcap.ff.oct$tot.SH), by=list(Mcap.ff.oct$vis), FUN=mean)
exp(means$x)  # S/H ratios in bleached vs. unbleached colonies
exp(means$x[1]) / exp(means$x[2])
# total vs. propD
plot(log(tot.SH) ~ asin(sqrt(propD)), data=Mcap.ff.oct,
     pch=c(24,21)[vis], bg=c("white", "black")[vis])
# • Figure 2: Relationship between symbiont community and bleaching -------------------------------
pdf("Figure2.pdf", height=3, width=3)
par(mfrow=c(1,2), mar=c(4,2,1,1), mgp=c(1.75,0.5,0))
# Plot barplot of C and D dominance in bleached and healthy corals
bars <- barplot(t(as.matrix(vis.chi$observed / rowSums(vis.chi$observed))), beside=F, xaxt="n", yaxt="n",
                ylab="", xlab="", cex.axis=0.75, cex.names=0.75, mgp=c(1, 0.25, 0), tck=-0.05)
axis(side=2, cex.axis=0.5, tck=-0.05, mgp=c(0.25,0.25,0))
mtext(side=2, text = "Proportion of colonies", line=1.1, cex=0.75)
text(x=bars, y=par("usr")[4] - 0.03, pos=3, xpd=T, cex=0.5,
     labels=c("n=30", "n=30"))
text(x=bars + 0.4, y=par("usr")[3] - 0.05, labels=c("Bleached", "Healthy"), cex=0.9, srt=45, pos=2, xpd=T)
legend(par("usr")[2] * c(0.95,1.15), c(0.8, 1.0), 
       legend=c("D", "C"), fill=c("gray95", "gray20"), xpd=NA, bty="n", cex=0.8, x.intersp=0.25)
# Plot bleaching severity in October by tdom:vis
par(mgp=c(2,0.5,0), mar=c(4,3.75,1,1))
plot(eff$fit, ylim=c(-5,-1.5), ylab="", yaxs="i", cex.axis=0.5,
     pch=21, bg="gray20", cex=2, line=1, bty="n", xpd=T, xaxt="n", xlab="", tck=-0.05, mgp=c(0.25,0.25,0))
mtext(side=2, "ln S/H", line=2, cex=0.75)
arrows(c(1,2,3), eff$fit - eff$se, c(1,2,3), eff$fit + eff$se, code=3, angle=90, length=0.075, xpd=T)
points(c(1,2,3), eff$fit, pch=21, bg=c("gray20", "gray20", "gray95"), cex=2, xpd=T)
axis(side=1, at=c(1,2,3), labels=NA, tck=-0.05)
text(c(1.5,2.5,3.5), par("usr")[3] - 0.3, xpd=T, srt=45, pos=2, cex=0.9,
     labels=c("Bleached (C)", "Healthy (C)", "Healthy (D)"))
dev.off()
# • Analysis: Symbiont community structure in each colony over time -------------------------------
clades <- melt(Mcap.f, id.vars=c("colony", "date", "vis", "reef", "tdom"), measure.vars="syms",
               factorsAsStrings=FALSE)
# Create matrix for image function
clades$value <- as.numeric(factor(clades$value))
clades <- dcast(clades, vis + colony + reef + tdom ~ date, drop=T)
clades[is.na(clades)] <- -1  # Recode missing values as -1
clades.m0 <- clades[with(clades, order(rev(vis), tdom, clades[, 5], clades[, 6], clades[, 7], 
                                    clades[, 8], clades[, 9], clades[, 10])), ]
clades.m <- as.matrix(clades.m0[,5:10], row.names=as.character(clades.m0$colony))
rownames(clades.m) <- as.character(clades$colony)
# How many colonies showed variable dominant symbionts over time
doms <- aggregate(Mcap.f$dom, by=list(colony=Mcap.f$colony), FUN=paste)
rownames(doms) <- doms$colony
domswitch <- lapply(doms$x, function(x) any(diff(as.numeric(factor(x)))!=0))  # TRUE if dom changed
domswitch <- data.frame(colony=doms$colony, domswitch=unlist(domswitch))
table(domswitch$domswitch)
df <- merge(domswitch, Mcap.f.summ)
model <- glm(domswitch ~ vis + reef + tdom, data=df, family=binomial)
dropterm(model, test="Chisq")
chisq.test(df$tdom, df$domswitch)$observed
chisq.test(df$vis, df$domswitch)$observed
# • Figure 3: Symbiont community structure in each colony over time ----------------------------------------
pdf("Figure3.pdf", width=3.5, height=6)
par(mfrow=c(1,1), mar=c(3,5,2,2), bg="white")
image(x=seq(1, ncol(clades.m)), y=seq(1, nrow(clades.m)), z=t(clades.m), 
      xaxt="n", yaxt="n", xlab="", ylab="",
      breaks=c(-1,0,1,2,3,4,5),
      col=c("white", rev(brewer.pal(11, "RdYlBu")[c(2,1,3,9,11)])))
# Plot date axis
axis(side=3, at=seq(1:6), labels=FALSE, cex.axis=0.75, par("tck"=-0.025), xpd=T)
text(1:6, par("usr")[4], xpd=T, cex=0.6, pos=3,
     labels=c("24 Oct\n2014", "4 Nov\n2014", "24 Nov\n2014", "16 Dec\n2014", "14 Jan\n2015", "6 May\n 2015"))
# Plot Bleached vs. Not Bleached rectangles
rect(par("usr")[1] - 2.25, par("usr")[3], par("usr")[1] - 1.25, par("usr")[4], xpd=T)
text(par("usr")[1] - 1.75, quantile(par("usr")[3:4])[c(2, 4)], labels=c("Not Bleached", "Bleached"),
     srt=90, xpd=2)
# Plot Row Side Colors
reefcols <- c("#bebada", "#8dd3c7", "#d9d9d9")
for (i in 1:nrow(clades)) {
  reef <- clades$reef[i]
  rect(par("usr")[1] - 1.25, par("usr")[3] + 1 * (i - 1), 
       par("usr")[1] - 0.25, par("usr")[3] + 1 * (i - 1) + 1, col=reefcols[as.numeric(reef)],
       xpd=T, border=NA)
}
rect(par("usr")[1] - 1.25, par("usr")[3], par("usr")[1] - 0.25, par("usr")[4], xpd=T)
lines(x=c(par("usr")[1] - 2.25, par("usr")[1] - 0.25), y=rep(quantile(par("usr")[3:4], 0.5), 2), xpd=T)
# Plot tdom side boxes
breaks <- c(0, which(diff(as.numeric(clades.m0$tdom))!=0), length(clades.m0$tdom))
doms <- clades.m0$tdom[breaks]
for (i in 2:length(breaks)) {
  rect(par("usr")[2] + 0.25, par("usr")[3] + breaks[i-1], par("usr")[2] + 0.75, par("usr")[3] + breaks[i], xpd=T)
}
for (i in 1:(length(breaks)-1)) {
  text(par("usr")[2] + 0.5, (breaks[i] + breaks[i+1]) / 2, paste(doms[i], "dominant"), xpd=T, srt=90,
       cex=0.75)
}
# Plot Row Side Color Key
for (i in 1:3) {
  rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1),
       par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1) - 1, xpd=T,
       border=NA, col=reefcols[i])
}
rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05, 
     par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - 3, xpd=T)
axis(side=2, xpd=T, pos=par("usr")[1] - 1.25, lwd=0, lwd.ticks=0,
     at=c(-1, -2, -3), labels=c("Rf 25", "Rf 44", "HIMB"), las=2, cex.axis=0.6, mgp=c(0,0.4,0))
# Plot Heatmap Key
for (i in 1:5) {
  rect(quantile(par("usr")[1:2], 0.1666 * i), quantile(par("usr")[3:4], 0) * -1.05,
       quantile(par("usr")[1:2], 0.1666 * (i + 1)), quantile(par("usr")[3:4], 0) * -1.05 - 1, xpd=T,
       border=NA, col=c(brewer.pal(11, "RdYlBu")[c(1,3,9,11)], "white")[i])
}
rect(quantile(par("usr")[1:2], 0.1666), quantile(par("usr")[3:4], 0) * -1.05, 
     quantile(par("usr")[1:2], 0.1666 * 6), quantile(par("usr")[3:4], 0) * -1.05 - 1, xpd=T)
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 0.75, pos=1, cex=0.6,
     x=seq(par("usr")[1], par("usr")[2], length=7)[-c(1,7)] + 0.5, 
     labels=c("D only", "D > C", "C > D", "C only", "no data"))
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 2.5, pos=1, cex=0.9,
     x=quantile(par("usr")[1:2], 0.5833),
     labels=expression(italic(Symbiodinium)~clades))
dev.off()
# • Analysis: Presence of background clade D and mixtures -----------------------------------------
# Presence of background D in bleached C colonies only over time
Cbleach <- Mcap.f[which(Mcap.f$tdom=="C" & Mcap.f$vis=="bleached"), ]
mod <- glmer(propD!=0 ~ fdate + (1|reef/colony), data=Cbleach, family=binomial)
mod
summary(mod)
# Look at presence of background clade D in bleached vs. non-bleached C colonies across all time points
Ccol <- Mcap.f[which(Mcap.f$tdom=="C"), ]
modr <- glmer(propD!=0 ~ vis * fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr2 <- glmer(propD!=0 ~ vis + fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr3 <- glmer(propD!=0 ~ fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr4 <- glmer(propD!=0 ~ (1|reef/colony), data=Ccol, family=binomial(link="logit"))
anova(modr, modr2, modr3, modr4)  # No significant effects of vis or date on presence of background clade D
# Look at frequency of CD mixtures vs. non-mixed frequency on diff dates
Mcap.f$mixed <- nchar(as.character(Mcap.f$syms)) - 1  # 1=mixed, 0=single clade
modr <- glmer(mixed ~ vis * fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr2 <- glmer(mixed ~ vis + fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr3 <- glmer(mixed ~ fdate + (1|reef/colony), data=Ccol, family=binomial(link="logit"))
modr4 <- glmer(mixed ~ (1|reef/colony), data=Ccol, family=binomial(link="logit"))
anova(modr, modr2, modr3, modr4)  # No significant effects of vis or date on presence of C/D mixtures
# =================================================================================================
# RECOVERY
# =================================================================================================
# • Analysis: Recovery dynamics --------------------------------------------------
#   Build piecewise polynomial model with knot at 82 days (January time point)
#   From October to January, fit a quadratic polynomial (1st element of degree=2)
#   From January to May, fit a linear model (2nd element of degree=1)
#   Function is continuous at time=82 days (smooth=0)
offset <- 0  # optional to center "days" axis at any point
sp <- function(x) gsp(x, knots=c(82 - offset), degree=c(2,1), smooth=0)
# Build and select model by backwards selection
mod.all.full <- lmerTest::lmer(log(Mcap.ff$tot.SH) ~ sp(Mcap.ff$days) * Mcap.ff$vis * Mcap.ff$tdom * Mcap.ff$reef + (1 | Mcap.ff$colony))
mod.all.full
modselect <- step(mod.all.full, lsmeans.calc=F, difflsmeans.calc=F)
# Print table of results
result <- modselect$anova.table
rownames(result) <- gsub("sp(Mcap.ff$days)", "time", fixed=T, rownames(result))
rownames(result) <- gsub("Mcap.ff$tdom", "clade", fixed=T, rownames(result))
rownames(result) <- gsub("Mcap.ff$reef", "reef", fixed=T, rownames(result))
rownames(result) <- gsub("Mcap.ff$vis", "bleach", fixed=T, rownames(result))
write.table(result, file="test.txt")
# Rebuild finalized model
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=Mcap.ff)
# Remove outliers with residuals > 2.5 s.d.'s from 0
out <- abs(residuals(mod.all)) > sd(residuals(mod.all)) * 2.5
Mcap.ff[out, ]  # outlying data points
# Refit model without outliers
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=Mcap.ff[!out, ])
summary(mod.all)  # look at parameter values
anova(mod.all)
pr <- summary(mod.all)$coefficients
write.csv(summary(mod.all)$coefficients, file="coef.csv")
# Generate predictions and confidence intervals using bootMer
pred.all <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")),
                         vis=factor(c("bleached", "not bleached")))
bootfit <- bootMer(mod.all, FUN=function(x) predict(x, pred.all, re.form=NA), nsim=999)
# Extract 95% confidence interval on predicted values
pred.all$fit <- predict(mod.all, pred.all, re.form=NA)
pred.all$lci <- apply(bootfit$t, 2, quantile, 0.025)
pred.all$uci <- apply(bootfit$t, 2, quantile, 0.975)
# Calculate when bleached are no longer different from healthy at each reef
data.frame(Effect(c("reef", "vis", "days"), mod.all, xlevels=list(days=c(0,82))))
plot(Effect(c("reef", "vis", "days"), mod.all, xlevels=list(days=1:194)), multiline=T, ci.style="bands")
plot(Effect(c("reef", "vis", "days"), mod.all2, xlevels=list(days=-81:112)), multiline=T, ci.style="bands")
summary(mod.all)
# ANALYZE MODEL COEFFICIENTS
summary(mod.all)
fixef(mod.all)
# SLOPES AT DAYS=0 -- test if each group is diff than zero
contrmat <- matrix(ncol=24, byrow=T, data=
                     c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  # 25.bleached == 0
                       0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  # 25.notbleached == 0
                       0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,  # 44.bleached == 0
                       0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,  # 44.notbleached == 0
                       0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,  # HIMB.bleached == 0
                       0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0)  # HIMB.notbleached == 0
                   )
rownames(contrmat) <- c("25.bleached","25.notbleached","44.bleached","44.notbleached",
                        "HIMB.bleached","HIMB.notbleached")
cust <- matrix(contrmat[1,] - contrmat[2,], nrow=1)
test <- glht(mod.all, linfct=contrmat)
test <- glht(mod.all, linfct=cust)
summary(test)

# Get model matrix defining each reef-vis group
mm <- unique(model.matrix(mod.all)[which(model.matrix(mod.all)[,"sp(days)D1(0)"]==82),])
rownames(mm) <- interaction(model.frame(mod.all)[rownames(mm), c("reef", "vis")])
mm[, grep("D1", colnames(mm), invert = T)] <- 0
mm <- (mm!=0) * 1
grp <- function(group) subset(mm, grepl(group, rownames(mm)))

test <- glht(mod.all, linfct=grp("HIMB.bleached")-grp("HIMB.not bleached"))
summary(test)
test <- glht(mod.all, linfct=grp("44.bleached")-grp("HIMB.not bleached"))
summary(test)
test <- glht(mod.all, linfct=mm)
summary(test)

slopecoefs <- names(fixef(mod.all))[grep("D1", names(fixef(mod.all)))]
slopecoefs[unique(c(1, grep("44", slopecoefs), grep("visnot bleached", slopecoefs)))]
c("44", "bleached") %in% slopecoefs
rownames(slopecoefs) <- rownames(summary(mod.all)$coefficients)[]
slopecoefs
names(fixef(mod.all))
null <- matrix(0, 1, length(fixef(mod.all)))

test <- glht(mod.all, linfct=null)
summary(test)

summary(mod.all)$coefficients[grep("D1", rownames(summary(mod.all)$coefficients)), ]
# HIMB not bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0), 1))
summary(test)
# HIMB bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# 44 not bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0), 1))
summary(test)
# 44 bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# 25 not bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# 25 bleached
test <- glht(mod.all, linfct=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)



# SLOPES AT DAY=150
newdat <- Mcap.ff[!out, ]
offset <- 150
newdat$days <- newdat$days - offset
mod.all2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=newdat)
#timeslopes <- summary(mod.all2)$coefficients[grep("D1", rownames(summary(mod.all2)$coefficients)), ]
summary(mod.all2)
### SLOPE OF POST_RECOV (JAN_MAY)
# Test 25 bleached == 0
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# Test 25 not bleached == 0
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# Test 25 bleached vs. 25 not bleached
test <- glht(mod.all, linfct=matrix(c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# Test 44 bleached == 0 ?????
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# Test 44 not bleached ==0
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0), 1))
summary(test)
# Test HIMB bleached == 0
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0), 1))
summary(test)
# Test HIMB not bleached == 0
test <- glht(mod.all, linfct=matrix(c(0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1), 1))
summary(test)




library(multcomp)
# test for diff between bleached and not bleached @ HIMB
test <- glht(mod.all2, linfct=matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,0), 1))
summary(test)
# test for diff between bleached @ 44 and bleached @ HIMB
test <- glht(mod.all2, linfct=matrix(c(0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)
# test for diff between bleached @ 44 and bleached @ HIMB
test <- glht(mod.all2, linfct=matrix(c(0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0), 1))
summary(test)






# # CAN MODEL BE REFIT SO THAT QUADRATIC EFFECT IS ONLY ALLOWED FOR BLEACHED CORALS?
# mf <- model.frame(mod.all)
# obj <- mf[, "sp(days)"]
# obj[which(mf$vis=="not bleached"), 2] <- 0
# mf$`sp(days)` <- obj
# mf
# df <- Mcap.ff[!out, ]
# mod.all2 <- lmerTest::lmer(log(df$tot.SH) ~ obj * df$vis * df$reef + (1 | df$colony))
# summary(mod.all2)
# model.frame(mod.all2)
# plot(effect(c("reef"), mod.all2))
# mod.all <- mod.all2
# • Figure 4: Recovery dynamics --------------------------------------------------
# Plotting function
plotreefs <- function(mod, pred) {
  dat <- model.frame(mod)
  dat <- merge(dat, unique(Mcap.ff[,c("colony", "reef", "vis")]))
  dat <- droplevels(dat)
  datsumm <- data.frame(expand.grid(reef=levels(dat$reef), vis=levels(dat$vis),
                                    days=as.numeric(as.character(levels(factor(dat$`sp(days)`[,1]))))),
                        mean=aggregate(dat$`log(tot.SH)`, by=list(interaction(dat$reef, dat$vis, dat$`sp(days)`[,1])), FUN=mean)$x,
                        sd=aggregate(dat$`log(tot.SH)`, by=list(interaction(dat$reef, dat$vis, dat$`sp(days)`[,1])), FUN=sd)$x,
                        se=aggregate(dat$`log(tot.SH)`, by=list(interaction(dat$reef, dat$vis, dat$`sp(days)`[,1])), 
                                     FUN=function(x) sd(x)/sqrt(length(x)))$x,
                        conf95=aggregate(dat$`log(tot.SH)`, by=list(interaction(dat$reef, dat$vis, dat$`sp(days)`[,1])), 
                                       FUN=function(x) sd(x)/sqrt(length(x)) * qt(0.975, length(x)-1))$x)
  datlist <- split(datsumm, f=datsumm$reef)
  datlist <- lapply(datlist, function(x) rev(split(x, f=x$vis)))
  predlist <- split(pred, f=pred$reef)
  predlist <- lapply(predlist, function(x) rev(split(x, f=x$vis)))
  par(mgp=c(1.75,0.4,0), oma=c(0,0,0,0))
  par(mar=c(0,3,0.3,1))
  for (reef in c("44", "25", "HIMB")) {
    with(datlist[[reef]], {
      # Create plot frame for each reef
      plot(NA, xlim=c(0,194), ylim=c(-6.5,-0.4), xaxt="n", bty="n", tck=-0.03, ylab="ln S/H")
      title(paste("Reef", reef), line=-0.9, adj=0, outer=F)
      # Plot model fit line and shaded CI for bleached and/or not bleached corals
      with(predlist[[reef]], {
        lapply(predlist[[reef]], function(vis) {
          addpoly(vis$days, vis$lci, vis$uci, col=alpha(reefcols[[reef]], 0.4), xpd=NA)
          lines(vis$days, vis$fit, lty=vislty[[vis$vis[1]]])
        })
      })
      # Plot raw data +/- standard deviation
      lapply(datlist[[reef]], function(vis) {
        arrows(vis$days, vis$mean + vis$se, vis$days, vis$mean - vis$se, code=3, angle=90, length=0.03, xpd=NA)
        points(vis$mean ~ vis$days, pch=vispch[[vis$vis[1]]], bg=visbg[[vis$vis[1]]], ylim=c(-7, 0.75))
      })
    })
    rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
  }
  axis(side=1, at=as.numeric(as.Date(c("2014-11-01", "2014-12-01", "2015-01-01", "2015-02-01", 
                                       "2015-03-01", "2015-04-01", "2015-05-01")) - as.Date("2014-10-24")),
       labels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"))
  return(list(predlist=predlist, datlist=datlist))
}
# Plot
#pdf("Figure4.pdf", width=3.5, height=7)
reefcols <- list(`25`="#bebada", `44`="#8dd3c7", HIMB="#d9d9d9")
vislty <- list("bleached"=2, "not bleached"=1)
vispch <- list("bleached"=24, "not bleached"=21)
visbg <- list("bleached"="white", "not bleached"="black")
layout(mat=matrix(c(1,2,3,4,4)))
modelplot <- plotreefs(mod.all, pred.all)
save1.x <- grconvertX(0, from='user', to='ndc' )
save1.y <- grconvertY(-6, from='user', to='ndc' )
save2.x <- grconvertX(82, from='user', to='ndc' )
save2.y <- grconvertY(-6, from='user', to='ndc' )
# Add zoom of bleached corals
par(mar=c(2,3,5,2))
plot(NA, ylim=c(-6,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n", ylab="", xlab="")
title("Recovery of\nbleached corals", adj=0, line=-2)
box(lty=3, col="black")
for (reef in c("44", "25", "HIMB")) {
  with(modelplot$predlist[[reef]]$bleached, {
    addpoly(days, lci, uci, col=alpha(reefcols[reef], 0.4))
    lines(days, fit, lty=2)
  })
}
# add zoom lines
segments(x0=0, y0=-1, x1=grconvertX(save1.x, from='ndc'), y1=grconvertY(save1.y, from='ndc'), lty=3, xpd=NA)
segments(x0=82, y0=-1, x1=grconvertX(save2.x, from='ndc'), y1=grconvertY(save2.y, from='ndc'), lty=3, xpd=NA)
dev.off()
# • Analysis: single time point statistics ----------------------
# October (bleached)
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
octmod <- lm(log(tot.SH) ~ reef * vis, data=Mcap.ff.oct)
anova(octmod)  # only vis is significant
octmod <- lm(log(tot.SH) ~ vis, data=Mcap.ff.oct)
octmeans <- data.frame(Effect("vis", octmod))$fit
anova(octmod)
# January (recovered)
Mcap.ff.jan <- Mcap.ff[which(Mcap.ff$fdate=="20150114"), ]
janmod <- lm(log(tot.SH) ~ reef * vis, data=Mcap.ff.jan)
summary(janmod)
anova(janmod)
janmean <- mean(log(Mcap.ff.jan$tot.SH))  # Grand geommean of january - no diff by reef or vis
# Comparisons
exp(octmeans[1])/exp(octmeans[2])  # Bleached are 88% lower than non-bleached in oct
exp(octmeans) / exp(janmean)  # Rel to jan, oct vals 95% lower in bleached, 61% lower in not-bleached
# Compare difference between OCt and Jan in each colony

# • Analysis: Recovery of neighbors of C vs. neighbors of D ----------------
C.nb <- as.numeric(as.character(unique(Mcap.ff.nb[which(Mcap.ff.nb$tdom=="C"), "sample"])))
D.nb <- as.numeric(as.character(unique(Mcap.ff.nb[which(Mcap.ff.nb$tdom=="D"), "sample"])))
Mcap.ff.b$neighbor <- factor(ifelse(Mcap.ff.b$sample %in% (C.nb - 1), "C", "D"))
Mcap.ff.b.octjan <- Mcap.ff.b[which(Mcap.ff.b$days %in% c(0, 11, 31, 52, 82)), ]
model <- lmerTest::lmer(log(tot.SH) ~ poly(days, 2) * neighbor * reef + (1 | sample), data=Mcap.ff.b.octjan)
summary(model)
dropterm(model, test="Chisq")
plot(Effect(c("days", "neighbor"), model), multiline=T)
plot(log(tot.SH) ~ days, col=reef, data=Mcap.ff.b.octjan)
# =================================================================================================
# • Temperature data ---------------------------------------------------------------
# Import temperature data
setwd("~/Documents/Academia/HIMB/K-Bay recovery 2014/Temp data/")
rf25.temp <- read.csv("Rf25_temps.csv")
rf44.temp <- read.csv("Rf44_temps.csv")
HIMB.temp <- read.csv("HIMB_temps.csv")
# Set date and time to POSIXct
rf25.temp$datetime <- paste(rf25.temp$date, rf25.temp$time)
rf25.temp$time <- as.POSIXct(rf25.temp$datetime, format="%m/%e/%y %H:%M:%S")
rf25.temp <- rf25.temp[, c(2,3,4)]
rf44.temp$datetime <- paste(rf44.temp$date, rf44.temp$time)
rf44.temp$time <- as.POSIXct(rf44.temp$datetime, format="%m/%e/%y %H:%M:%S")
rf44.temp <- rf44.temp[, c(2,3,4)]
HIMB.temp$datetime <- paste(HIMB.temp$date, HIMB.temp$time)
HIMB.temp$time <- as.POSIXct(HIMB.temp$datetime, format="%m/%e/%y %H:%M:%S")
HIMB.temp <- HIMB.temp[, c(2,3,4)]
# Subset data only up until threshdate
threshdate <- as.POSIXct("2015-09-01", format="%F")
rf25.temp <- rf25.temp[which(rf25.temp$time < threshdate), ]
rf44.temp <- rf44.temp[which(rf44.temp$time < threshdate), ]
HIMB.temp <- HIMB.temp[which(HIMB.temp$time < threshdate), ]

# Aggregate temperature data by daily mean, minimum, and maximum
rf25 <- aggregate(data.frame(mean=rf25.temp$temp_cal), by=list(date=as.Date(rf25.temp$time)), FUN=mean)
rf25$min <- aggregate(rf25.temp$temp_cal, by=list(as.Date(rf25.temp$time)), FUN=min)$x
rf25$max <- aggregate(rf25.temp$temp_cal, by=list(as.Date(rf25.temp$time)), FUN=max)$x
rf44 <- aggregate(data.frame(mean=rf44.temp$temp_cal), by=list(date=as.Date(rf44.temp$time)), FUN=mean)
rf44$min <- aggregate(rf44.temp$temp_cal, by=list(as.Date(rf44.temp$time)), FUN=min)$x
rf44$max <- aggregate(rf44.temp$temp_cal, by=list(as.Date(rf44.temp$time)), FUN=max)$x
HIMB <- aggregate(data.frame(mean=HIMB.temp$temp_cal), by=list(date=as.Date(HIMB.temp$time)), FUN=mean)
HIMB$min <- aggregate(HIMB.temp$temp_cal, by=list(as.Date(HIMB.temp$time)), FUN=min)$x
HIMB$max <- aggregate(HIMB.temp$temp_cal, by=list(as.Date(HIMB.temp$time)), FUN=max)$x

# Plot daily mean, min, and max for each reef
plot(mean ~ date, rf44, type="l", col="darkgreen", ylim=c(21.5,30))
polygon(x=c(rf44$date, rev(rf44$date)), y=c(rf44$min, rev(rf44$max)),
        border=NA, col=alpha("darkgreen", 0.2))
lines(mean ~ date, HIMB, type="l", col="red")
polygon(x=c(HIMB$date, rev(HIMB$date)), y=c(HIMB$min, rev(HIMB$max)),
        border=NA, col=alpha("red", 0.2))
rf25.split <- split(rf25, f=rf25$date < as.Date("2014-11-21", format="%F"))
lines(mean ~ date, rf25.split[[1]], type="l", col="blue")
polygon(x=c(rf25.split[[1]]$date, rev(rf25.split[[1]]$date)), y=c(rf25.split[[1]]$min, rev(rf25.split[[1]]$max)),
        border=NA, col=alpha("blue", 0.2))
lines(mean ~ date, rf25.split[[2]], type="l", col="blue")
polygon(x=c(rf25.split[[2]]$date, rev(rf25.split[[2]]$date)), y=c(rf25.split[[2]]$min, rev(rf25.split[[2]]$max)),
        border=NA, col=alpha("blue", 0.2))
legend("bottomright", lty=1, col=c("darkgreen", "blue", "red"), legend=c("44","25","HIMB"))
