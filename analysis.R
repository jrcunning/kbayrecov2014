# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 27 August, 2015

# Run setup script
source("setup.R")
# Set seed
set.seed(39059978)

# =================================================================================================
# DOMINANT SYMBIONTS AND BLEACHING
# =================================================================================================
# • Analysis: Symbiont clades in all samples and colonies across time points ----------------------
# Presence/absence of C and D
symtab <- table(Mcap.f$syms)
symtab2 <- c(symtab[1], symtab[2] + symtab[3], symtab[4])
par(mfrow=c(2,1), mar=c(3,5,2,1))
barplot(symtab2, ylab="Number of samples", xlab="Symbiont clades detected")
clades <- aggregate(Mcap.f$syms, by=list(sample=Mcap.f$sample), FUN=paste, collapse="")
clades$C[grep("C", clades$x)] <- "C"
clades$D[grep("D", clades$x)] <- "D"
clades$present <- ifelse(is.na(clades$C), ifelse(is.na(clades$D), "none", "D only"), ifelse(is.na(clades$D), "C only", "C+D"))
barplot(table(clades$present), ylab="Number of colonies", xlab="Symbiont clades detected")
compare <- data.frame(Colonies=matrix(prop.table(table(clades$present))), Samples=prop.table(symtab2))
# Proportion clade D
logDC <- Mcap.f$logDC
DC <- exp(logDC)
propD <- DC/(DC + 1)
propD[is.na(propD)] <- 1
# • Figure: overall symbiont community composition in Montipora capitata  ---------------
symcols <- rev(brewer.pal(11, "RdYlBu")[c(2,1,3,9,11)])
par(mfrow=c(1,2), mar=c(2,3,2,1), mgp=c(1.75,0.25,0), lwd=1)
b1 <- barplot(as.matrix(compare), beside=F, ylab="Proportion", space=0.5, 
              col = c("gray20", "gray", "white"), cex.names = 1, cex.axis=0.8)
text(b1, par("usr")[4], pos=3, xpd=T,
     labels = paste("n=", c(length(clades$present), sum(symtab2)), sep=""))
save1.x <- grconvertX(par("usr")[2], from='user', to='ndc')
save1.y <- grconvertY(compare[1,2], from='user', to='ndc')
save2.y <- grconvertY(compare[1,2] + compare[2,2], from="user", to="ndc")
legend(x=c(3.25,6), y=c(0,0.2), legend=c("D only", "C and D", "C only"), xpd=NA, cex=1,
       fill=c("white", "gray", "gray20"), bty="n")
par(mar=c(10,2,2,4))
hist(propD[which(propD > 0 & propD < 1)], breaks=c(0,0.25,0.5,0.75,1), plot=T,
     main="", xaxt="n", xlab="", yaxt="n", ylab="", col = "gray")
box()
axis(side=1, at=seq(0,1,0.25), labels=c(">0","","0.5","","<1"), cex.axis=0.75)
mtext(side=1, "Proportion clade D", line=1.5, cex=0.75)
axis(side=4, cex.axis=0.75)
mtext(side=4, "Number of samples", line=1.5, cex=0.75)
segments(x0 = par("usr")[1], y0 = par("usr")[3], 
         x1 = grconvertX(save1.x, from='ndc'), 
         y1 = grconvertY(save1.y, from='ndc'), lty=1, xpd=NA)
segments(x0 = par("usr")[1], y0 = par("usr")[4], 
         x1 = grconvertX(save1.x, from='ndc'), 
         y1 = grconvertY(save2.y, from='ndc'), lty=1, xpd=NA)
# • Analysis and Figure: Statistical tests of symbiont clade, reef, and visual status --------------------
# Summarize data - reef, visual status, and dominant symbiont for each colony
Mcap.f.summ <- unique(Mcap.f[, c("sample", "vis", "reef", "tdom")])
# Logistic regression testing effect of visual appearance and reef on dominant clade
model <- glm(tdom ~ vis * reef, data=Mcap.f.summ, family=binomial(link="logit"))
anova(model, test="Chisq")  # only vis is significant
model <- glm(tdom ~ vis, data=Mcap.f.summ, family=binomial(link="logit"))
anova(model, test="Chisq")
par(mfrow=c(1,2), mar=c(5,4,2,2), mgp=c(2.5,0.5,0))
plot(tdom ~ vis, data=Mcap.f.summ,  # 100% of bleached colonies were clade C; 57% of notbleached were D
     axes=F, ylab="Proportion of colonies", xlab="Visual diagnosis")
axis(side=2, line=0.3)
axis(side=1, at=quantile(par("usr")[1:2], c(0.25,0.75)), labels=c("Bleached", "Healthy"), cex.axis=1)
legend(c(1,1.2),c(0.8,1.0), legend=c("D", "C"), fill=c("gray90", "gray30"), xpd=NA, bty="n", cex=0.8, x.intersp=0.25)
#CONCLUSION: bleached are all clade C, not bleached may be clade C or D 50/50
# Bleaching severity in October -- effect of symbiont clade and visual appearance
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
model <- lm(log(tot.SH) ~ tdom:vis, data=Mcap.ff.oct)
anova(model, test="F")  # tdom:vis significant, (reef is not significant if included)
#plot(effect("tdom:vis", model))
eff <- data.frame(effect("tdom:vis", model))[c(1,3,4), ]
TukeyHSD(aov(model))
# Conclusion: Bleached C colonies had lowest S/H, Unbleached C and D colonies had same S/H
par(mgp=c(2,0.5,0), mar=c(5,5,2,3))
plot(eff$fit, ylim=c(-6,-1), ylab="", yaxs="i",
     pch=21, bg="gray20", cex=2, line=1, bty="n", xpd=T, xaxt="n", xlab="")
mtext(side=2, "ln S/H", line=3)
arrows(c(1,2,3), eff$lower, c(1,2,3), eff$upper, code=3, angle=90, length=0.1, xpd=T)
axis(side=1, at=c(1,2,3), labels=NA)
text(c(1,2,3), par("usr")[3] - 0.5, xpd=T, srt=45, pos=1,
     labels=c("Bleached (C)", "Healthy (C)", "Healthy (D)"))
# • Analysis: Symbiont clades in each colony over time ------------------------------------------------------
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
# • Figure: Symbiont clades in each colony over time ----------------------------------------
par(mfrow=c(1,1), mar=c(5,7,3,3), bg="white")
image(x=seq(1,ncol(clades.m)), y=seq(1,nrow(clades.m)), z=t(clades.m), 
      xaxt="n", yaxt="n", xlab="", ylab="",
      breaks=c(-1,0,1,2,3,4,5),
      col=c("white", rev(brewer.pal(11, "RdYlBu")[c(2,1,3,9,11)])))
# Plot date axis
axis(side=3, at=seq(1:6), labels=FALSE, cex.axis=0.75, par("tck"=-0.025), xpd=T)
text(1:6, par("usr")[4], xpd=T, cex=0.75, pos=3,
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
  text(par("usr")[2] + 0.5, (breaks[i] + breaks[i+1]) / 2, paste(doms[i], "dominant"), xpd=T, srt=90)
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
     at=c(-1, -2, -3), labels=c("Rf 25", "Rf 44", "HIMB"), las=2, cex.axis=0.75, mgp=c(0,0.4,0))
# Plot Heatmap Key
for (i in 1:5) {
  rect(quantile(par("usr")[1:2], 0.1666 * i), quantile(par("usr")[3:4], 0) * -1.05,
       quantile(par("usr")[1:2], 0.1666 * (i + 1)), quantile(par("usr")[3:4], 0) * -1.05 - 1, xpd=T,
       border=NA, col=c(brewer.pal(11, "RdYlBu")[c(1,3,9,11)], "white")[i])
}
rect(quantile(par("usr")[1:2], 0.1666), quantile(par("usr")[3:4], 0) * -1.05, 
     quantile(par("usr")[1:2], 0.1666 * 6), quantile(par("usr")[3:4], 0) * -1.05 - 1, xpd=T)
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 0.75, pos=1, cex=0.7,
     x=seq(par("usr")[1], par("usr")[2], length=7)[-c(1,7)] + 0.5, 
     labels=c("D only", "D > C", "C > D", "C only", "no data"))
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 2, pos=1, cex=0.9,
     x=quantile(par("usr")[1:2], 0.5833),
     labels=expression(italic(Symbiodinium)~clades))
# • Analysis: D:C ratio in each colony over time --------------------------------------------------
# Rearrange data: colonies as rows, logDC ratio for each date as columns, reef + vis as columns
logDC <- melt(Mcap.f, id.vars=c("sample", "date", "vis", "reef"), measure.vars=c("logDC"),
                 factorsAsStrings=FALSE)
logDC <- dcast(logDC, reef + vis + sample ~ date, drop=T)
logDC <- logDC[with(logDC, order(vis, reef, sample)), ]
rownames(logDC) <- logDC$sample
# Order colonies (rows) by bleached/notbleached, then amount of D at each date
logDC <- logDC[with(logDC, order(rev(vis), logDC[,4], logDC[,5], logDC[,6], logDC[,7],
                                 logDC[,8], logDC[,9])), ]
logDC.m <- as.matrix(logDC[,4:9])
rownames(logDC.m) <- rownames(logDC)
# Find range of non-infinite ratios
range(logDC.m[which(!is.infinite(logDC.m))], na.rm=T)
# Assign infinite values to just beyond this range (-Inf=-12, +Inf=+12)
logDC.m[which(logDC.m==-Inf)] <- -12
logDC.m[which(logDC.m==Inf)] <- 12
# • Figure: D:C ratio in each colony over time -----------------------------------------------------
par(mfrow=c(1,1), mar=c(5,7,3,1), bg="white")
image(x=seq(1,ncol(logDC.m)), y=seq(1,nrow(logDC.m)), z=t(logDC.m), 
      xaxt="n", yaxt="n", xlab="", ylab="",
      breaks=c(-12, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12),
      col=rev(brewer.pal(10, "RdYlBu")))
# Plot date axis
axis(side=3, at=seq(1:6), labels=FALSE, cex.axis=0.75, par("tck"=-0.025), xpd=T)
text(1:6, par("usr")[4], xpd=T, cex=0.75, pos=3,
     labels=c("24 Oct\n2014", "4 Nov\n2014", "24 Nov\n2014", "16 Dec\n2014", "14 Jan\n2015", "6 May\n 2015"))
# Plot Bleached vs. Not Bleached rectangles
rect(par("usr")[1] - 2.25, par("usr")[3], par("usr")[1] - 1.25, par("usr")[4], xpd=T)
text(par("usr")[1] - 1.75, quantile(par("usr")[3:4])[c(2, 4)], labels=c("Not Bleached", "Bleached"),
     srt=90, xpd=2)

# Plot Row Side Colors
reefcols <- c("#bebada", "#8dd3c7", "#d9d9d9")
for (i in 1:nrow(logDC.m)) {
  reef <- logDC$reef[i]
  rect(par("usr")[1] - 1.25, par("usr")[3] + 1 * (i - 1), 
       par("usr")[1] - 0.25, par("usr")[3] + 1 * (i - 1) + 1, col=reefcols[as.numeric(reef)],
       xpd=T, border=NA)
}
rect(par("usr")[1] - 1.25, par("usr")[3], par("usr")[1] - 0.25, par("usr")[4], xpd=T)
lines(x=c(par("usr")[1] - 2.25, par("usr")[1] - 0.25), y=rep(quantile(par("usr")[3:4], 0.5), 2), xpd=T)

# Plot Row Side Color Key
for (i in 1:3) {
  rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1),
       par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1) - 1, xpd=T,
       border=NA, col=reefcols[i])
}
rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05, 
     par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - 3, xpd=T)
axis(side=2, xpd=T, pos=par("usr")[1] - 1.25, lwd=0, lwd.ticks=0,
     at=c(-1, -2, -3), labels=c("Rf 25", "Rf 44", "HIMB"), las=2, cex.axis=0.75, mgp=c(0,0.4,0))

# Plot heatmap Key
breaks=c(-12, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12)
keycol=c(rev(brewer.pal(10, "RdYlBu")), "white")
for (i in 1:(length(breaks))) {
  rect(par("usr")[1] + (i-1) * diff(par("usr")[1:2]) / (length(breaks)), quantile(par("usr")[3:4], 0) * -1.05, 
       par("usr")[1] + i * diff(par("usr")[1:2]) / (length(breaks)), quantile(par("usr")[3:4], 0) * -1.05 - 1,
       col= keycol[i], xpd=T, border = NA)
}
rect(par("usr")[1], quantile(par("usr")[3:4], 0) * -1.05, par("usr")[2], quantile(par("usr")[3:4], 0) * -1.05 - 1, xpd=T)

axis(side=1, xpd=T, pos=quantile(par("usr")[3:4], 0) * -1.05 - 1,
     at=seq(par("usr")[1], par("usr")[2], length=length(breaks) + 1)[2:10],
     labels=breaks[2:10], cex.axis=0.5, mgp=c(0,-0.2,0), tck=-0.01)

text(x=quantile(par("usr")[1:2], 0.5), y=quantile(par("usr")[3:4], 0) * -1.05 - 1.5, 
     labels = "log D:C", cex=0.75, xpd=T, pos=1)
lines(x=rep(quantile(par("usr")[1:2], 0.045), 2), 
      y=c(quantile(par("usr")[3:4], 0) * -1.05 - 0.5, quantile(par("usr")[3:4], 0) * -1.05 - 2.5), 
      xpd=T)
lines(x=rep(quantile(par("usr")[1:2], 0.9545), 2), 
      y=c(quantile(par("usr")[3:4], 0) * -1.05 - 0.5, quantile(par("usr")[3:4], 0) * -1.05 - 2.5), 
      xpd=T)
lines(x=rep(quantile(par("usr")[1:2], 0.8636), 2), 
      y=c(quantile(par("usr")[3:4], 0) * -1.05 - 0.5, quantile(par("usr")[3:4], 0) * -1.05 - 2.5), 
      xpd=T)
text(x=c(quantile(par("usr")[1:2], c(0.045, 0.8636, 0.9545))),
     y=quantile(par("usr")[3:4], 0) * -1.05 - 3.25,
     labels = c("C only", "D only", "no data"), xpd=T, cex=0.5)

# -------------------------------------------------------------------------------------------------
# =================================================================================================
# RECOVERY
# =================================================================================================
# • Bleached corals
# -------------------------------------------------------------------------------------------------
# - Figure: How bleached were corals in October-----------
mod <- lm(log(tot.SH) ~ tdom:vis, data=Mcap.ff)
dropterm(mod, test="F")
plot(effect("tdom:vis", mod))
# - Figure: symbiont abundance over time by reef (plot raw data) ----------------
Mcap.ff.b <- subset(Mcap.ff, vis=="bleached")
# Visualize raw data / trajectories for each colony over time
xyplot(log(tot.SH) ~ days | reef, groups = ~ sample, data = Mcap.ff.b[order(Mcap.ff.b$days), ],
       type = "o", layout=c(3, 1), main="bleached colonies")
# - Analysis: symbiont abundance over recovery period only (October to January) -------------------
# Subset data: bleached corals between october and january
Mcap.ff.b <- subset(Mcap.ff, vis=="bleached")
Mcap.ff.b.octjan <- subset(Mcap.ff.b, date <= "2015-01-14")
# Build models with linear, quadratic, and cubic functions of time
lmodel <- lme4::lmer(log(tot.SH) ~ days * reef + (days|sample), data=Mcap.ff.b.octjan)
qmodel <- lme4::lmer(log(tot.SH) ~ poly(days, 2) * reef + (poly(days,2)|sample), data=Mcap.ff.b.octjan)
cmodel <- lme4::lmer(log(tot.SH) ~ poly(days, 3) * reef + (poly(days,3)|sample), data=Mcap.ff.b.octjan)
anova(lmodel, qmodel, cmodel)  # quadratic model is best
# Build quadratic model with random intercept only
qmodel2 <- lme4::lmer(log(tot.SH) ~ poly(days, 2) * reef + (1|sample), data=Mcap.ff.b.octjan)
anova(qmodel, qmodel2)  # model with random intercept only is best
mod.b.octjan <- qmodel2  # select model with random intercept
dropterm(mod.b.octjan, test="Chisq")  # LRT indicates interaction term is significant (compared to model without interaction)
summary(mod.b.octjan)  # Model parameters
mcp.fnc(mod.b.octjan)  # Model diagnostics -- residuals are homoscedastic, QQplot is normal
# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(mod.b.octjan, Mcap.ff.b.octjan, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
mod.b.octjan <- update(mod.b.octjan, data = rm.outliers$data)
# Generate model predictions and confidence intervals using bootMer
newdat.b.octjan <- expand.grid(days=seq(0,82,1), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(mod.b.octjan, FUN=function(x) predict(x, newdat.b.octjan, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat.b.octjan$pred <- predict(mod.b.octjan, newdat.b.octjan, re.form=NA)
newdat.b.octjan$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat.b.octjan$uci <- apply(bootfit$t, 2, quantile, 0.95)
# - Figure: symbiont abundance over recovery period only -------------------------------------------
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(2,2,2,2))
with(newdat.b.octjan[newdat.b.octjan$reef=="HIMB", ], {
  plot(pred ~ days, type="l", col="red", ylim=c(-6,-1))
  addpoly(days, lci, uci, col=alpha("red", 0.5))
})
with(newdat.b.octjan[newdat.b.octjan$reef=="25", ], {
  lines(pred ~ days, col="blue")
  addpoly(days, lci, uci, col=alpha("blue", 0.5))
})
with(newdat.b.octjan[newdat.b.octjan$reef=="44", ], {
  lines(pred ~ days, col="darkgreen")
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.5))
})
# - Analysis: symbiont abundance over full time series (October to May) ---------------------------
# Subset data
Mcap.ff.b <- subset(Mcap.ff, vis=="bleached")
# Build piecewise polynomial model with knot at 82 days (January time point)
#   From October to January, fit a quadratic polynomial (1st element of degree=2)
#   From January to May, fit a linear model (2nd element of degree=1)
#   Function is continuous at time=82 days (smooth=0)
sp <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
# Build models with random intercept and slope, and random intercept only
mod.b <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (1 | sample), data=Mcap.ff.b)
mod.b2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (sp(days) | sample), data=Mcap.ff.b)
anova(mod.b, mod.b2)  # Model with random intercept only is best
# Test parameters
summary(mod.b)  # Coefficients for sp(days)D1(0) indicate slope at days=0
                # Coefficients for sp(days)D2(0) are neg. for convex, pos. for concave
lstrends(mod.b, ~ reef, var = "days")  # Where do these values come from?
dropterm(mod.b, test="Chisq")
# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(mod.b, Mcap.ff.b, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
mod.b <- update(mod.b, data = rm.outliers$data)
# Generate predictions and confidence intervals using bootMer
newdat.b <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(mod.b, FUN=function(x) predict(x, newdat.b, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat.b$pred <- predict(mod.b, newdat.b, re.form=NA)
newdat.b$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat.b$uci <- apply(bootfit$t, 2, quantile, 0.95)
newdat.b <- split(newdat.b, f=newdat.b$reef)
# Summarize raw data for plotting
mdf <- model.frame(mod.b)
model.data.summ <- data.frame(expand.grid(reef=levels(mdf$reef), days=as.numeric(as.character(levels(factor(mdf$`sp(days)`[,1]))))),
                              mean=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=sd)$x)
model.data.summ <- split(model.data.summ, f=model.data.summ$reef)
# - Figure: symbiont abundance over full time series (mean ± SD) with model fit and C.I.'s --------
layout(mat=matrix(c(1,2,3,4,4)))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1), mar=c(0,2,0,0))
with(model.data.summ$`44`, {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.b$"44", lines(days, pred))
  with(newdat.b$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`25`, {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.b$"25", lines(days, pred))
  with(newdat.b$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`HIMB`, {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-7,1), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.b$"HIMB", lines(days, pred))
  with(newdat.b$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-6,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat.b$"HIMB", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.b$"25", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.b$"44", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# -------------------------------------------------------------------------------------------------
# • Non-bleached corals
# -------------------------------------------------------------------------------------------------
# - Figure: symbiont abundance in non-bleached corals over time ---------------------------------
# Subset data
Mcap.ff.nb <- subset(Mcap.ff, vis=="not bleached")
Mcap.ff.nb <- droplevels(subset(Mcap.ff.nb, tdom!="CD"))  # Remove colony 40
xyplot(log(tot.SH) ~ days | reef + tdom, groups = ~ sample, data=Mcap.ff.nb[order(Mcap.ff.nb$days), ],
       type = "o", layout=c(3, 2), main="not bleached colonies")
# - Analysis: symbiont abundance over full time series (October to May) --------------------
# Subset data
Mcap.ff.nb <- subset(Mcap.ff, vis=="not bleached")
Mcap.ff.nb <- droplevels(subset(Mcap.ff.nb, tdom!="CD"))  # Remove mixed colony 40 for now
# Build model
spq <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
mod.nb <- lmerTest::lmer(log(tot.SH) ~ spq(days) * reef + (1 | sample), data=Mcap.ff.nb)
mod.nb2 <- lmerTest::lmer(log(tot.SH) ~ spq(days) * reef + (spq(days) | sample), data=Mcap.ff.nb)
anova(mod.nb, mod.nb2)  # Model with random intercept only is best
dropterm(mod.nb, test="Chisq")
sp <- function(x) gsp(x, knots=c(82), degree=c(1,1), smooth=0)
mod.nb3 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (1 | sample), data=Mcap.ff.nb)
anova(mod.nb, mod.nb3) # Quadratic fit is not better, p = 0.1979
mod.nb4 <- lmerTest::lmer(log(tot.SH) ~ sp(days) + reef + (1 | sample), data=Mcap.ff.nb)
anova(mod.nb3, mod.nb4)  # model with reef*days interaction is almost sig. by LRT and has lower AIC--> keep
mod.nb <- mod.nb3  # Choose linear model with random intercept only

# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(mod.nb, Mcap.ff.nb, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
mod.nb <- update(mod.nb, data = rm.outliers$data)
# Generate predictions and prediction intervals using bootMer over days 0-194
newdat.nb <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(mod.nb, FUN=function(x) predict(x, newdat.nb, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat.nb$pred <- predict(mod.nb, newdat.nb, re.form=NA)
newdat.nb$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat.nb$uci <- apply(bootfit$t, 2, quantile, 0.95)
#newdat.nb <- split(newdat.nb, f=interaction(newdat.nb$reef, newdat.nb$tdom))
newdat.nb <- split(newdat.nb, f=newdat.nb$reef)
# Plot data
mdf.nb <- model.frame(mod.nb)
mdf.nb.summ <- data.frame(expand.grid(reef=levels(mdf.nb$reef),
                                          days=as.numeric(as.character(levels(factor(mdf.nb$`sp(days)`[,1]))))),
                              mean=aggregate(mdf.nb$`log(tot.SH)`, by=list(interaction(mdf.nb$reef, mdf.nb$`sp(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf.nb$`log(tot.SH)`, by=list(interaction(mdf.nb$reef, mdf.nb$`sp(days)`[,1])), FUN=sd)$x)

mdf.nb.summ <- split(mdf.nb.summ, f=mdf.nb.summ$reef)
# - Figure: symbiont abundance over full time series --------------
layout(mat=matrix(c(1,2,3,4,4)))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1))
par(mar=c(0,2,0,0))
with(mdf.nb.summ$"44", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.nb$"44", lines(days, pred))
  with(newdat.nb$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -5, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(mdf.nb.summ$"25", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"25", lines(days, pred))
  with(newdat.nb$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -5, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(mdf.nb.summ$"HIMB", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-7,1), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"HIMB", lines(days, pred))
  with(newdat.nb$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -5, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-5,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat.nb$"HIMB", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.nb$"25", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.nb$"44", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# # Plot D's
# par(mar=c(0,2,0,0))
# with(model.data.summ$"44.D", {
#   plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
#   arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
#   with(newdat.nb$"44.D", lines(days, pred))
#   with(newdat.nb$"44.D", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
#   rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
# })
# par(mar=c(0,2,0,0))
# with(model.data.summ$"25.D", {
#   plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
#   arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
#   with(newdat.nb$"25.D", lines(days, pred))
#   with(newdat.nb$"25.D", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
#   rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
# })
# par(mar=c(0,2,0,0))
# with(model.data.summ$"HIMB.D", {
#   plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
#   arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
#   with(newdat.nb$"HIMB.D", lines(days, pred))
#   with(newdat.nb$"HIMB.D", addpoly(days, lci, uci, col=alpha("red", 0.3)))
#   rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
# })
# par(mar=c(1,2,5,0))
# plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
# mtext(side=3, text = "Days", line=2.5, cex=0.75)
# box(lty=2, col=alpha("black", 0.8))
# with(newdat.nb$"HIMB.D", {
#   lines(days, pred)
#   addpoly(days, lci, uci, col=alpha("red", 0.3))
# })
# with(newdat.nb$"25.D", {
#   lines(days, pred)
#   addpoly(days, lci, uci, col=alpha("blue", 0.3))
# })
# with(newdat.nb$"44.D", {
#   lines(days, pred)
#   addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
# })







# -------------------------------------------------------------------------------------------------
# • All corals
# -------------------------------------------------------------------------------------------------
# Build piecewise polynomial model with knot at 82 days (January time point) -----------
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
# PLOT all corals----------------------------------------
layout(mat=matrix(c(1,2,3,4,4,5,6,7,8,8), ncol=2))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1), mar=c(0,2,0,0))
with(model.data.summ$"44.bleached", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.bleached", lines(days, pred))
  with(newdat.all$"44.bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25.bleached", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.bleached", lines(days, pred))
  with(newdat.all$"25.bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB.bleached", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-7,1), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.bleached", lines(days, pred))
  with(newdat.all$"HIMB.bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-6,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat.all$"HIMB.bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.all$"25.bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.all$"44.bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# Plot not bleached
par(mar=c(0,2,0,0))
with(model.data.summ$"44.not bleached", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.not bleached", lines(days, pred))
  with(newdat.all$"44.not bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25.not bleached", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-7,1), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.not bleached", lines(days, pred))
  with(newdat.all$"25.not bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB.not bleached", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-7,1), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.not bleached", lines(days, pred))
  with(newdat.all$"HIMB.not bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-6,-1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat.all$"HIMB.not bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.all$"25.not bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.all$"44.not bleached", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# Figure: Plot bleached and non-bleached on same panels------------
pdf("recovplot2.pdf", width=3.5, height=7)
layout(mat=matrix(c(1,2,3,4,4)))
par(mgp=c(1.75,0.4,0), oma=c(0,0,0,0))
par(mar=c(0,3,0,1))
with(model.data.summ$"44.bleached", {
  plot(mean ~ days, pch=21, bg=alpha("darkgreen", 0.2), ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03, 
       ylab="ln S/H")
  title("Reef 44", line=-1.5, adj=0.9)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.bleached", lines(days, pred, lty=2))
  with(newdat.all$"44.bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"44.not bleached", {
  points(mean ~ days, pch=21, bg="darkgreen", ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.not bleached", lines(days, pred))
  with(newdat.all$"44.not bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
})
par(mar=c(0,3,0,1))
with(model.data.summ$"25.bleached", {
  plot(mean ~ days, pch=21, bg=alpha("blue", 0.2), ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03, ylab="ln S/H")
  title("Reef 25", line=-1.5, adj=0.9)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.bleached", lines(days, pred, lty=2))
  with(newdat.all$"25.bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"25.not bleached", {
  points(mean ~ days, pch=21, bg="blue", ylim=c(-7,0.75), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.not bleached", lines(days, pred))
  with(newdat.all$"25.not bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
})
par(mar=c(0,3,0,1))
with(model.data.summ$"HIMB.bleached", {
  plot(mean ~ days, pch=21, bg=alpha("red", 0.2), ylim=c(-7,0.75), bty="n", tck=-0.03, xaxt="n", ylab="ln S/H")
  title("HIMB", line=-1.5, adj=0.9)
  axis(side=1, at=as.numeric(as.Date(c("2014-11-01", "2014-12-01", "2015-01-01", "2015-02-01", 
                                       "2015-03-01", "2015-04-01", "2015-05-01")) - as.Date("2014-10-24")),
       labels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"))
  
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.bleached", lines(days, pred, lty=2))
  with(newdat.all$"HIMB.bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
})
with(model.data.summ$"HIMB.not bleached", {
  points(mean ~ days, pch=21, bg="red", ylim=c(-7,0.75), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.not bleached", lines(days, pred))
  with(newdat.all$"HIMB.not bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
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
  lines(days, pred, lty=2)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.all$"25.bleached", {
  lines(days, pred, lty=2)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.all$"44.bleached", {
  lines(days, pred, lty=2)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# add zoom lines
segments(x0=0, y0=-1, x1=grconvertX(save1.x, from='ndc'), y1=grconvertY(save1.y, from='ndc'), lty=3, xpd=NA)
segments(x0=82, y0=-1, x1=grconvertX(save2.x, from='ndc'), y1=grconvertY(save2.y, from='ndc'), lty=3, xpd=NA)
dev.off()














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
