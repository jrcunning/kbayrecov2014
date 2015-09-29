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

# Analyze patterns in clade composition across reefs
cladecomp <- merge(unique(Mcap.f[, c("colony", "reef", "tdom")]), colonies)
# Analyze proportions of mixed colonies across reefs
chisq1 <- chisq.test(cladecomp$reef, cladecomp$present)  # no differences
chisq1$observed; chisq1$p.value  # No differences
# Analyze proportions of dominant symbionts across reefs
chisq2 <- chisq.test(cladecomp$reef, cladecomp$tdom)
chisq2$observed; chisq2$p.value  # No differences
# • Figure 1: Overall symbiont community composition in Montipora capitata  ---------------
#pdf(file="Figure1.pdf", height=3, width=3)
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
#dev.off()
# • Analysis: Relationship between dominant symbiont clade and bleaching -------------------------------
# Summarize data - reef, visual status, and dominant symbiont for each colony
Mcap.f.summ <- unique(Mcap.f[, c("colony", "vis", "reef", "tdom")])
# Logistic regression testing effect of visual appearance and reef on dominant clade
model <- glm(tdom ~ vis * reef, data=Mcap.f.summ, family=binomial(link="logit"))
anova(model, test="Chisq")  # only vis is significant
# • Figure 2: Relationship between symbiont community and bleaching -------------------------------
#pdf("Figure2.pdf", height=3, width=3)
par(mfrow=c(1,2), mar=c(4,3,1,1), mgp=c(1.75,0.5,0))
plot(tdom ~ vis, data=Mcap.f.summ,  # 100% of bleached colonies were clade C; 57% of notbleached were D
     axes=F, ylab="Proportion of colonies", xlab="", cex.axis=0.75, mgp=c(0.25,0.25,0),
     col=c("gray20", "gray95"))
axis(side=2, line=0.3, cex.axis=0.75, tck=-0.05, mgp=c(0.5,0.25,0))
text(x=quantile(par("usr")[1:2], c(0.25, 0.75)), y=par("usr")[4] - 0.03, pos=3, xpd=T, cex=0.75,
     labels=c("n=30", "n=30"))
text(x=quantile(par("usr")[1:2], c(0.4,0.9)), y=par("usr")[3] - 0.05, 
     labels=c("Bleached", "Healthy"), cex.axis=1, srt=45, pos=2, xpd=T)
legend(c(1,1.2),c(0.8,1.0), legend=c("D", "C"), fill=c("gray95", "gray20"), xpd=NA, bty="n", cex=0.8, x.intersp=0.25)
#CONCLUSION: bleached are all clade C, not bleached may be clade C or D 50/50
# Bleaching severity in October -- effect of symbiont clade and visual appearance
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
model <- lm(log(tot.SH) ~ tdom:vis, data=Mcap.ff.oct)
anova(model, test="F")  # tdom:vis significant, (reef is not significant if included)
#plot(effect("tdom:vis", model))
eff <- data.frame(effect("tdom:vis", model))[c(1,3,4), ]
TukeyHSD(aov(model))
# Conclusion: Bleached C colonies had lowest S/H, Unbleached C and D colonies had same S/H
par(mgp=c(2,0.5,0), mar=c(4,4,1,1))
plot(eff$fit, ylim=c(-6,-1), ylab="", yaxs="i", cex.axis=0.75,
     pch=21, bg="gray20", cex=2, line=1, bty="n", xpd=T, xaxt="n", xlab="", tck=-0.05, mgp=c(0.5,0.25,0))
mtext(side=2, "ln S/H", line=2.5)
arrows(c(1,2,3), eff$lower, c(1,2,3), eff$upper, code=3, angle=90, length=0.1, xpd=T)
points(c(1,2,3), eff$fit, pch=21, bg=c("gray20", "gray20", "gray95"), cex=2, xpd=T)
axis(side=1, at=c(1,2,3), labels=NA, tck=-0.05)
text(c(1.5,2.5,3.5), par("usr")[3] - 0.3, xpd=T, srt=45, pos=2, cex=0.9,
     labels=c("Bleached (C)", "Healthy (C)", "Healthy (D)"))
#dev.off()
# • Analysis: Symbiont community structure in each colony over time ------------------------------------------------------
clades <- melt(Mcap.f, id.vars=c("sample", "date", "vis", "reef", "tdom"), measure.vars="syms",
               factorsAsStrings=FALSE)
# Create matrix for image function
clades$value <- as.numeric(factor(clades$value))
clades <- dcast(clades, vis + sample + reef + tdom ~ date, drop=T)
clades[is.na(clades)] <- -1  # Recode missing values as -1
clades <- clades[with(clades, order(rev(vis), tdom, clades[, 5], clades[, 6], clades[, 7], 
                                    clades[, 8], clades[, 9], clades[, 10])), ]
clades.m <- as.matrix(clades[,5:10], row.names=as.character(clades$sample))
rownames(clades.m) <- as.character(clades$sample)
# • Figure 3: Symbiont community structure in each colony over time ----------------------------------------
pdf("Figure3.pdf", width=3.5, height=6)
par(mfrow=c(1,1), mar=c(3,5,2,2), bg="white")
image(x=seq(1,ncol(clades.m)), y=seq(1,nrow(clades.m)), z=t(clades.m), 
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
clades$tdom
breaks <- c(0, which(diff(as.numeric(clades$tdom))!=0), length(clades$tdom))
doms <- clades$tdom[breaks]
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
# Look at presence of background clade D in bleached vs. non-bleached C colonies across all time points
Ccol <- Mcap.f[which(Mcap.f$tdom=="C"), ]
modr <- glmer(propD!=0 ~ vis * fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr2 <- glmer(propD!=0 ~ vis + fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr3 <- glmer(propD!=0 ~ fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr4 <- glmer(propD!=0 ~ (1|reef/sample), data=Ccol, family=binomial(link="logit"))
anova(modr, modr2, modr3, modr4)  # No significant effects of vis or date on presence of background clade D
# Look at frequency of CD mixtures vs. non-mixed frequency on diff dates
Mcap.f$mixed <- nchar(as.character(Mcap.f$syms)) - 1  # 1=mixed, 0=single clade
modr <- glmer(mixed ~ vis * fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr2 <- glmer(mixed ~ vis + fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr3 <- glmer(mixed ~ fdate + (1|reef/sample), data=Ccol, family=binomial(link="logit"))
modr4 <- glmer(mixed ~ (1|reef/sample), data=Ccol, family=binomial(link="logit"))
anova(modr, modr2, modr3, modr4)  # No significant effects of vis or date on presence of C/D mixtures
# =================================================================================================
# RECOVERY
# =================================================================================================
# • Analysis: Recovery dynamics at different reefs --------------------------------------------------
#   Build piecewise polynomial model with knot at 82 days (January time point)
#   From October to January, fit a quadratic polynomial (1st element of degree=2)
#   From January to May, fit a linear model (2nd element of degree=1)
#   Function is continuous at time=82 days (smooth=0)
sp <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
# Build models with random intercept and slope, and random intercept only
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | sample), data=Mcap.ff)
mod.all2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (sp(days) | sample), data=Mcap.ff)
anova(mod.all, mod.all2)  # Model with random intercept only is best
# Test significance of factors in model
dropterm(mod.all, test="Chisq")
# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(mod.all, Mcap.ff, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
mod.all <- update(mod.all, data = rm.outliers$data)
# Generate predictions and confidence intervals using bootMer
newdat.all <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")),
                          vis=factor(c("bleached", "not bleached")))
bootfit <- bootMer(mod.all, FUN=function(x) predict(x, newdat.all, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat.all$pred <- predict(mod.all, newdat.all, re.form=NA)
newdat.all$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat.all$uci <- apply(bootfit$t, 2, quantile, 0.95)
newdat.all <- split(newdat.all, f=interaction(newdat.all$reef, newdat.all$vis))
# Summarize raw data for plotting
mdf <- model.frame(mod.all)
model.data.summ <- data.frame(expand.grid(reef=levels(mdf$reef), vis=levels(mdf$vis),
                                          days=as.numeric(as.character(levels(factor(mdf$`sp(days)`[,1]))))),
                              mean=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$vis, mdf$`sp(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$vis, mdf$`sp(days)`[,1])), FUN=sd)$x)
model.data.summ <- split(model.data.summ, f=interaction(model.data.summ$reef, model.data.summ$vis))
# • Figure 4: Recovery dynamics at different reefs --------------------------------------------------
reefcols <- c("#bebada", "#8dd3c7", "#d9d9d9")
reefcols <- c("blue", "darkgreen", "red")
pdf("recovplot2.pdf", width=3.5, height=7)
layout(mat=matrix(c(1,2,3,4,4)))
par(mgp=c(1.75,0.4,0), oma=c(0,0,0,0))
par(mar=c(0,3,0,1))
with(model.data.summ$"44.bleached", {
  plot(mean ~ days, type="n", ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03, 
       ylab="ln S/H")
  title("Reef 44", line=-1.5, adj=0.9)
  with(newdat.all$"44.bleached", addpoly(days, lci, uci, col=alpha(reefcols[2], 0.3)))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.bleached", lines(days, pred, lty=2))
  points(mean ~ days, pch=24, bg=reefcols[2], ylim=c(-7,0.75), bty="n", tck=-0.03)
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"44.not bleached", {
  with(newdat.all$"44.not bleached", addpoly(days, lci, uci, col=alpha(reefcols[2], 0.3)))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.not bleached", lines(days, pred))
  points(mean ~ days, pch=21, bg=reefcols[2], ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03)
})
par(mar=c(0,3,0,1))
with(model.data.summ$"25.bleached", {
  plot(mean ~ days, type="n", ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03, ylab="ln S/H")
  title("Reef 25", line=-1.5, adj=0.9)
  with(newdat.all$"25.bleached", addpoly(days, lci, uci, col=alpha(reefcols[1], 0.3)))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.bleached", lines(days, pred, lty=2))
  points(mean ~ days, pch=24, bg=reefcols[1], ylim=c(-7,0.75), bty="n", tck=-0.03)
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"25.not bleached", {
  with(newdat.all$"25.not bleached", addpoly(days, lci, uci, col=alpha(reefcols[1], 0.3)))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.not bleached", lines(days, pred))
  points(mean ~ days, pch=21, bg=reefcols[1], ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03)
})
par(mar=c(0,3,0,1))
with(model.data.summ$"HIMB.bleached", {
  plot(mean ~ days, type="n", ylim=c(-7,0.75), bty="n", tck=-0.03, xaxt="n", ylab="ln S/H")
  title("HIMB", line=-1.5, adj=0.9)
  axis(side=1, at=as.numeric(as.Date(c("2014-11-01", "2014-12-01", "2015-01-01", "2015-02-01", 
                                       "2015-03-01", "2015-04-01", "2015-05-01")) - as.Date("2014-10-24")),
       labels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"))
  
  with(newdat.all$"HIMB.bleached", addpoly(days, lci, uci, col=alpha(reefcols[3], 0.3)))
  with(newdat.all$"HIMB.bleached", lines(days, pred, lty=2))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  points(mean ~ days, pch=24, bg=reefcols[3], ylim=c(-7,0.75), bty="n", tck=-0.03)
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"HIMB.not bleached", {
  with(newdat.all$"HIMB.not bleached", addpoly(days, lci, uci, col=alpha(reefcols[3], 0.3)))
  with(newdat.all$"HIMB.not bleached", lines(days, pred))
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  points(mean ~ days, pch=21, bg=reefcols[3], ylim=c(-7,0.75), bty="n", tck=-0.03)
})
save1.x <- grconvertX(0, from='user', to='ndc' )
save1.y <- grconvertY(-6, from='user', to='ndc' )
save2.x <- grconvertX(82, from='user', to='ndc' )
save2.y <- grconvertY(-6, from='user', to='ndc' )
par(mar=c(2,3,5,2))
plot(NA, ylim=c(-6,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n", ylab="", xlab="")
title("Close-up: recovery of\nbleached corals", adj=0, line=-2)
#mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=3, col="black")
with(newdat.all$"HIMB.bleached", {
  addpoly(days, lci, uci, col=alpha(reefcols[3], 0.3))
  lines(days, pred, lty=2)
})
with(newdat.all$"25.bleached", {
  addpoly(days, lci, uci, col=alpha(reefcols[1], 0.3))
  lines(days, pred, lty=2)
})
with(newdat.all$"44.bleached", {
  addpoly(days, lci, uci, col=alpha(reefcols[2], 0.3))
  lines(days, pred, lty=2)
})
# add zoom lines
segments(x0=0, y0=-1, x1=grconvertX(save1.x, from='ndc'), y1=grconvertY(save1.y, from='ndc'), lty=3, xpd=NA)
segments(x0=82, y0=-1, x1=grconvertX(save2.x, from='ndc'), y1=grconvertY(save2.y, from='ndc'), lty=3, xpd=NA)
dev.off()
# • Analysis: Recovery of neighbors of C vs. neighbors of D
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
