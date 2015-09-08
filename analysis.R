# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 27 August, 2015

# Run setup script
source("setup.R")
# Set seed
set.seed(39059978)

# =================================================================================================
# PATTERNS IN DOMINANT SYMBIONTS AND BLEACHING
# =================================================================================================
# • Analysis: Symbiont clades in each colony over time ------------------------------------------------------
clades <- melt(Mcap.f, id.vars=c("sample", "date", "vis", "reef"), measure.vars="syms",
               factorsAsStrings=FALSE)
# Create matrix for image function
clades$value <- as.numeric(factor(clades$value))
clades <- dcast(clades, vis + sample + reef ~ date, drop=T)
clades[is.na(clades)] <- -1  # Recode missing values as -1
clades <- clades[with(clades, order(rev(vis), clades[, 4], clades[, 5], clades[, 6], 
                                    clades[, 7], clades[, 8], clades[, 9])), ]
clades.m <- as.matrix(clades[,4:9], row.names=as.character(clades$sample))
rownames(clades.m) <- as.character(clades$sample)
# • Figure: Symbiont clades in each colony over time ----------------------------------------
par(mar=c(5,7,3,1), bg="white")
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
# • Statistical tests of symbiont clade, reef, and visual status --------------------
# NEED TO CORRECT FOR COPY NUMBER BEFORE ANY OF THIS IS VALID
# Summarize data - reef, visual status, and dominant symbiont for each colony
Mcap.f.oct <- Mcap.f[which(Mcap.f$fdate=="20141024"), ]
Mcap.f.summ <- unique(Mcap.f.oct[, c("sample", "vis", "reef", "dom", "tdom")])

# Figure: dominant symbiont clades of bleached and unbleached colonies at each reef
par(mfrow=c(1, 3), mar=c(3,3,3,3))    #USE TDOM OR Oct-DOM???????  NEED TO CORRECT FOR COPY # FIRST.
plot(dom ~ vis, Mcap.f.summ[which(Mcap.f.summ$reef=="HIMB"), ],
     xlab="visual", ylab="dominant clade", main="HIMB")
plot(dom ~ vis, Mcap.f.summ[which(Mcap.f.summ$reef=="25"), ],
     xlab="visual", ylab="", main="25")
plot(dom ~ vis, Mcap.f.summ[which(Mcap.f.summ$reef=="44"), ],
     xlab="visual", ylab="", main="44")

# Logistic regression testing effect of visual appearance and reef on dominant clade
model <- glm(dom ~ vis + reef, data=Mcap.f.summ, family=binomial(link="logit"))
anova(model, test="Chisq")  # vis and reef are significant, but not interaction
par(mfrow=c(1,2))
plot(tdom ~ vis, data=Mcap.f.summ)  # 100% of bleached colonies were clade C; 57% of notbleached were D
plot(tdom ~ reef, data=Mcap.f.summ)  # HIMB has more clade D colonies (9/20) than 25&44 (4/20 each)
#CONCLUSION: bleached are all clade C, not bleached may be clade C or D
#            reefs 25 and 44 80% likely to be clade C, at HIMB, nearly 58/42

# Follow up: look at only not-bleached colonies
# Ho: Not-bleached colonies' dominant symbiont did not depend on reef
# Ha: The dominant symbiont in not-bleached colonies depended on reef
notbleached <- Mcap.f.summ[which(Mcap.f.summ$vis=="not bleached"), ]
par(mfrow=c(1,1))
plot(tdom ~ reef, data=notbleached, main="Not-bleached colonies")
model <- glm(tdom ~ reef, data=notbleached, family=binomial(link="logit"))
dropterm(model, test="Chisq")
# Conclusion: dominant symbiont in non-bleached colonies differed among reefs
#             Non-bleached at HIMB were 90% D, Non-bleached at 25 and 44 were 60% C
          # Interpretation: stress was not as severe at 25 and 44, allowing many C colonies to remain unbleached
          #                 whereas at HIMB, only D's avoided bleaching
          #             OR  host-mediated resistance more common at 25 and 44


##Ho: Visual appearance did not depend on symbiont clade
##Ha: Visual appearance was linked to symbiont clade 
par(mfrow=c(1,3))
plot(vis ~ tdom, Mcap.f.summ[which(Mcap.f.summ$reef=="HIMB"), ],
     xlab="dominant clade", ylab="visual", main="HIMB")
plot(vis ~ tdom, Mcap.f.summ[which(Mcap.f.summ$reef=="25"), ],
     xlab="dominant clade", ylab="", main="25")
plot(vis ~ tdom, Mcap.f.summ[which(Mcap.f.summ$reef=="44"), ],
     xlab="dominant clade", ylab="", main="44")
# Logistic regression testing effect of clade and reef on visual appearance
model <- glm(vis ~ tdom * reef, data=Mcap.f.summ, family=binomial(link="logit"))
anova(model, test="Chisq")  # symbiont clade is significant
par(mfrow=c(1,1))
plot(vis ~ tdom, data=Mcap.f.summ)
# CONCLUSION: visual bleaching depends on clade--70% of clade C colonies appeared bleached, while
#             no clade D colonies appeared bleached
# -------------------------------------------------------------------------------------------------

# =================================================================================================
# QUANTITATIVE ANALYSIS OF SYMBIONT ABUNDANCE
# =================================================================================================
# • Bleached corals
# -------------------------------------------------------------------------------------------------
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
  plot(pred ~ days, type="l", col="red", ylim=c(-3.25,1.5))
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
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.b$"44", lines(days, pred))
  with(newdat.b$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`25`, {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.b$"25", lines(days, pred))
  with(newdat.b$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`HIMB`, {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.b$"HIMB", lines(days, pred))
  with(newdat.b$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
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

# - Analysis: Piecewise linear regression over full time period (October to May) --------------------
# Subset data
Mcap.ff.nb <- subset(Mcap.ff, vis=="not bleached")
Mcap.ff.nb <- droplevels(subset(Mcap.ff.nb, tdom!="CD"))  # Remove mixed colony 40 for now
# Build model
sp <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
mod.nb <- lmerTest::lmer(log(tot.SH) ~ sp(days) * tdom * reef + (1 | sample), data=Mcap.ff.nb)
mod.nb2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * tdom * reef + (sp(days) | sample), data=Mcap.ff.nb)
anova(mod.nb, mod.nb2)  # Model with random intercept only is best
splin <- function(x) gsp(x, knots=c(82), degree=c(1,1), smooth=0)
mod.nb3 <- lmerTest::lmer(log(tot.SH) ~ splin(days) * tdom * reef + (1 | sample), data=Mcap.ff.nb)
anova(mod.nb, mod.nb3)  # Quadratic fit is better
mod.nb <- mod.nb  # Choose quadratic model with random intercept only
dropterm(mod.nb, test="Chisq")  # FULL INT NOT SIG
mod.nb4 <- lmerTest::lmer(log(tot.SH) ~ sp(days) + tdom + reef + 
                            sp(days):tdom + sp(days):reef + tdom:reef +
                            (1 | sample), data=Mcap.ff.nb)
mod.nb5 <- lmerTest::lmer(log(tot.SH) ~ sp(days) + tdom + reef + 
                            sp(days):tdom + sp(days):reef +
                            (1 | sample), data=Mcap.ff.nb)
mod.nb6 <- lmerTest::lmer(log(tot.SH) ~ sp(days) + reef + 
                            (1 | sample), data=Mcap.ff.nb)
mod.nb7 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (1 | sample), data=Mcap.ff.nb)
anova(mod.nb6, mod.nb7)
mod.nb8 <- lmerTest::lmer(log(tot.SH) ~ splin(days) + reef + (1 | sample), data=Mcap.ff.nb)
anova(mod.nb6, mod.nb8)
mod.nb <- mod.nb8
dropterm(mod.nb, test="Chisq")
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
mdf <- model.frame(mod.nb)
model.data.summ <- data.frame(expand.grid(reef=levels(mdf$reef),
                                          days=as.numeric(as.character(levels(factor(mdf$`sp(days)`[,1]))))),
                              mean=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=sd)$x)

model.data.summ <- split(model.data.summ, f=model.data.summ$reef)
# PLOTTT
layout(mat=matrix(c(1,2,3,4,4)))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1))
par(mar=c(0,2,0,0))
with(model.data.summ$"44", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.nb$"44", lines(days, pred))
  with(newdat.nb$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"25", lines(days, pred))
  with(newdat.nb$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"HIMB", lines(days, pred))
  with(newdat.nb$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
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
# Plot D's
par(mar=c(0,2,0,0))
with(model.data.summ$"44.D", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.nb$"44.D", lines(days, pred))
  with(newdat.nb$"44.D", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25.D", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"25.D", lines(days, pred))
  with(newdat.nb$"25.D", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB.D", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.nb$"HIMB.D", lines(days, pred))
  with(newdat.nb$"HIMB.D", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat.nb$"HIMB.D", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat.nb$"25.D", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat.nb$"44.D", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})







# -------------------------------------------------------------------------------------------------
# • All corals
# -------------------------------------------------------------------------------------------------
# Build piecewise polynomial model with knot at 82 days (January time point)
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
#PLOT
layout(mat=matrix(c(1,2,3,4,4,5,6,7,8,8), ncol=2))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1), mar=c(0,2,0,0))
with(model.data.summ$"44.bleached", {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.bleached", lines(days, pred))
  with(newdat.all$"44.bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25.bleached", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.bleached", lines(days, pred))
  with(newdat.all$"25.bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB.bleached", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.bleached", lines(days, pred))
  with(newdat.all$"HIMB.bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
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
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat.all$"44.not bleached", lines(days, pred))
  with(newdat.all$"44.not bleached", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"25.not bleached", {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"25.not bleached", lines(days, pred))
  with(newdat.all$"25.not bleached", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$"HIMB.not bleached", {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat.all$"HIMB.not bleached", lines(days, pred))
  with(newdat.all$"HIMB.not bleached", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
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










# Bleaching severity in October -- effect of symbiont clade and visual appearance------------
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
model <- lm(log(tot.SH) ~ dom:vis, data=Mcap.ff.oct)
anova(model, test="F")  # tdom:vis significant, (reef is not significant if included)
plot(effect("dom:vis", model))
TukeyHSD(aov(model))  # BLEACHED C SAME AS NONBLEACHED D--- BECAUSE OF COPY NUMBER?
# Conclusion: Bleached C colonies had lowest S/H, Unbleached C and D colonies had same S/H







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
