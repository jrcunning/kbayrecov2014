# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 27 August, 2015
#
# Run setup script
source("setup.R")
# Set seed
set.seed(39059978)

# =================================================================================================
# SYMBIONT COMMUNITY STRUCTURE AND BLEACHING
# =================================================================================================
# • Analysis: Symbiodinium community structure --------------------------
# Proportion of samples with C only, D only, and C+D mixtures
symtab <- table(Mcap.f$syms)
samples <- c(symtab[1], symtab[2] + symtab[3], symtab[4])
prop.table(samples)
# Proportion clade D in samples with C+D mixtures
propD <- Mcap.f$propD[which(Mcap.f$propD > 0 & Mcap.f$propD < 1)]
hist(propD)
range(propD)
# Percent of samples with >10% non-dominant symbiont (between 10% and 90% clade D)
sum(prop.table(hist(propD, plot=F)$counts)[2:9])
# Proportion of colonies with C only, D only, and C+D mixtures, aggregated over time
colonies <- aggregate(Mcap.f$syms, by=list(colony=Mcap.f$colony), FUN=paste, collapse="")
colonies$C[grep("C", colonies$x)] <- "C"
colonies$D[grep("D", colonies$x)] <- "D"
colonies$present <- ifelse(is.na(colonies$C), ifelse(is.na(colonies$D), "none", "D only"), ifelse(is.na(colonies$D), "C only", "C+D"))
prop.table(table(colonies$present))
# Summarize clade composition of samples and colonies
clades <- data.frame(Colonies=matrix(prop.table(table(colonies$present))), Samples=prop.table(samples))
clades
# Proportion of colonies with overall C or D dominance (most abundant over time)
propdom <- prop.table(table(Mcap.f[which(Mcap.f$fdate=="20141024"), "tdom"]))
propdom
# • Figure 2: Occurrence of Symbiodinium clades C and D in M. capitata samples and colonies  ---------------
pdf(file="Figure2*.pdf", height=3, width=3)
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

# • Analysis: Influence of Symbiodinium on bleaching response -------------------------------
# Test for differences in dominant clade between bleached and unbleached colonies
Mcap.f.summ <- unique(Mcap.f[, c("colony", "vis", "reef", "tdom")])
vis.chi <- chisq.test(Mcap.f.summ$vis, Mcap.f.summ$tdom)
vis.chi[c("p.value", "observed")]
# Test for differences in S/H between C-bleached, C-notbleached, and D-notbleached
sh.mod <- lm(log(tot.SH) ~ tdom:vis, data=Mcap.ff.oct)
sh.lsm <- lsmeans(sh.mod, specs=c("tdom", "vis"))[c(1,3,4)] 
cld(sh.lsm)
bvnb <- summary(lsmeans(update(sh.mod, ~ vis), specs="vis")) # compare bleached vs. not-bleached
exp(bvnb$lsmean)  # S/H in bleached vs. non-bleached corals
# Compare presence of D in C-dominated colonies that bleached vs. did not bleach
Ccol <- subset(Mcap.f, tdom=="C")
model <- glmer(!is.na(D.SH) ~ vis + (1|colony), family="binomial", data=Ccol)
lsmeans(model, specs = "vis", contr="pairwise")
# Test for differences in S/H among reefs in bleached corals
Mcap.ff.oct.b <- subset(Mcap.ff, fdate=="20141024" & vis=="bleached")
bleach <- lm(log(tot.SH) ~ reef, data=Mcap.ff.oct.b)
anova(bleach, test="F")
# Test for differences in dominant clades of unbleached colonies across reefs
with(Mcap.f.summ[which(Mcap.f.summ$vis=="not bleached"), ], {
  chisq.test(reef, tdom)[c("observed", "p.value")]
})  # No difference
# • Figure 3: Influence of Symbiodinium on bleaching in M. capitata -------------------------------
pdf("Figure3*.pdf", height=3.18898, width=3.18898)  # 81 mm square
par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1.75,0.5,0))
# Plot barplot of C and D dominance in bleached and healthy corals
bars <- barplot(t(as.matrix(vis.chi$observed / rowSums(vis.chi$observed))), beside=F, xaxt="n", yaxt="n", yaxs="i",
                ylab="", xlab="", cex.axis=0.75, cex.names=0.75, mgp=c(1, 0.25, 0), tck=-0.05)
axis(side=2, cex.axis=0.5, tck=-0.05, mgp=c(0.25,0.25,0))
mtext(side=2, text = "Proportion of colonies", line=1.1, cex=0.75)
text(x=bars, y=par("usr")[4] - 0.03, pos=3, xpd=T, cex=0.5,
     labels=c("n=30", "n=30"))
text(x=bars, y=par("usr")[3], labels=c("B", "NB"), cex=0.75, pos=1, xpd=T)
legend(par("usr")[2] * c(0.95,1.15), c(0.5, 0.9), 
       legend=c("D", "C"), fill=c("gray95", "gray20"), xpd=NA, bty="n", cex=0.7, x.intersp=0.25)
text(x=-0.25, y=par("usr")[4], labels="A", xpd=T, pos=2)
# Plot bleaching severity in October by tdom:vis
par(mgp=c(2,0.5,0), mar=c(2,3.5,1,1))
dat <- summary(sh.lsm)
plot(dat$lsmean, ylim=c(-5,-1.5), ylab="", yaxs="i", cex.axis=0.5,
     line=1, bty="n", xpd=T, xaxt="n", xlab="", tck=-0.05, mgp=c(0.25,0.25,0))
mtext(side=2, "ln S/H", line=2, cex=0.75)
arrows(c(1,2,3), dat$lsmean - dat$SE, c(1,2,3), dat$lsmean + dat$SE, code=3, angle=90, length=0.075, xpd=T)
points(c(1,2,3), dat$lsmean, pch=21, bg=c("gray20", "gray20", "gray95"), cex=2, xpd=T)
text(c(1,2,3), dat$lsmean + dat$SE, labels=c("a","b","b"), pos=3, xpd=T, cex=0.5)
axis(side=1, at=c(1,2,3), labels=NA, tck=0.05)
text(c(1,2,3), par("usr")[3], xpd=T, pos=1, cex=0.75,
     labels=c("B\n(C)", "NB\n(C)", "NB\n(D)"))
text(x=-0.125, y=par("usr")[4], labels="B", xpd=T, pos=2)
dev.off()

# • Analysis: Temporal paterns in Symbiodinium composition -------------------------------
# Determine occurrence of dominant symbiont stability vs. shuffling (at any point in time)
doms <- melt(Mcap.f, id.vars=c("colony", "date", "vis", "reef", "tdom"), measure.vars="dom",
             factorsAsStrings=FALSE)
doms <- dcast(doms, vis + colony + reef + tdom ~ date, drop=T)
domswitch <- apply(doms[,5:10], 1, function(x) any(diff(as.numeric(factor(x[!is.na(x)])))!=0)) 
domswitch <- data.frame(colony=doms$colony, domswitch=unlist(domswitch))
domswitch <- merge(doms, domswitch); domswitch   # TRUE if dominant symbiont changed at any time
table(domswitch$domswitch)  # Number of colonies that showed change in dominant symbiont
# Test whether shuffling (at any time point) is related to bleaching, reef, or dominant symbiont
model <- glm(domswitch ~ vis * reef * tdom, data=domswitch, family=binomial)
anova(model, test="Chisq")  # No significant effects of vis, reef, or tdom
# Test whether shuffling events are related to date and/or bleaching
shuff <- data.frame(t(apply(doms[,5:10], 1, FUN=function(x) diff(as.numeric(factor(x))))))  # CtoD=1, DtoC=-1
colnames(shuff) <- levels(Mcap.f$fdate)[-1]
rowSums(shuff)  # if 0, shuffling was transient within study period (4 colonies =/= 0)
shuff <- cbind(doms[,1:4], ifelse(shuff==0,0,1))  # 1 if shuffled, 0 if stable at each time point
shuff.m <- melt(shuff, id.vars=c("colony", "vis", "reef", "tdom"), value.name="shuff", variable.name="fdate")
shuff.mod <- glmer(shuff ~ vis * fdate + (1|colony), data=shuff.m, family="binomial")
dropterm(shuff.mod, test="Chisq")  # vis*fdate interaction not significant
shuff.mod <- update(shuff.mod, ~ vis + fdate + (1|colony))  # remove interaction term
dropterm(shuff.mod, test="Chisq")  # test vis and fdate main effects
# Test whether presence of background symbionts is related to date and/or bleaching
Ccol <- subset(Mcap.f, tdom=="C")
bgD <- glmer(is.na(D.SH) ~ vis * fdate + (1|colony), data=Ccol, family="binomial")
dropterm(bgD, test="Chisq")  # vis*fdate interaction not significant
bgD <- update(bgD, ~ vis + fdate + (1|colony))  # remove interaction term
dropterm(bgD, test="Chisq") # vis and fdate are not significant
Dcol <- subset(Mcap.f, tdom=="D")
bgC <- glmer(is.na(C.SH) ~ fdate + (1|colony), data=Dcol, family="binomial")
dropterm(bgC, test="Chisq")  # fdate not significant
# • Figure 4: Symbiont community structure in each colony over time ----------------------------------------
# Create matrix for image function
clades <- melt(Mcap.f, id.vars=c("colony", "date", "vis", "reef", "tdom"), measure.vars="syms",
               factorsAsStrings=FALSE)
clades$value <- as.numeric(factor(clades$value))
clades <- dcast(clades, vis + colony + reef + tdom ~ date, drop=T)
clades[is.na(clades)] <- -1  # Recode missing values as -1
clades[which(clades$colony=="129"), 8:10] <- -2  # Recode mortality as -2
clades.m0 <- clades[with(clades, order(rev(vis), tdom, clades[, 5], clades[, 6], clades[, 7], 
                                       clades[, 8], clades[, 9], clades[, 10])), ]
clades.m <- as.matrix(clades.m0[,5:10])
rownames(clades.m) <- as.character(clades.m0$colony)
# Plot figure
pdf("Figure4*.pdf", width=3.5, height=6)
par(mfrow=c(1,1), mar=c(3,5,2,2), bg="white")
image(x=seq(1, ncol(clades.m)), y=seq(1, nrow(clades.m)), z=t(clades.m), 
      xaxt="n", yaxt="n", xlab="", ylab="",
      breaks=c(-2,-1.1,0,1,2,3,4,5),
      col=c("black", "white", rev(brewer.pal(11, "RdYlBu")[c(2,1,3,9,11)])))
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
for (i in 1:nrow(clades.m0)) {
  reef <- clades.m0$reef[i]
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
x <- quantile(par("usr")[1:2], probs=seq(0, 1, length.out=7))
y <- rep(quantile(par("usr")[3:4], 0) * -1.05, 2) - c(0, 1)
rect(x[1], y[1], x[7], y[2], xpd=T)
for (i in 1:6) {
  rect(x[i], y[1], x[i + 1], y[2], xpd=T,
       border=NA, col=c(brewer.pal(11, "RdYlBu")[c(1,3,9,11)], "white", "black")[i])
}


text(xpd=T, y=y[1] - 0.75, pos=1, cex=0.6,
     x=seq(par("usr")[1], par("usr")[2], length=7)[-7] + 0.5, 
     labels=c("D only", "D > C", "C > D", "C only", "no data", "dead"))
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 2.5, pos=1, cex=0.9,
     x=quantile(par("usr")[1:2], 0.5),
     labels=expression(italic(Symbiodinium)~clades))
dev.off()



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
# Build full model with fixed effects of vis, tdom, reef, and time, random effect of colony
mod.all.full <- lmerTest::lmer(log(Mcap.ff$tot.SH) ~ sp(Mcap.ff$days) * Mcap.ff$vis * Mcap.ff$tdom * Mcap.ff$reef + (1 | Mcap.ff$colony))
# Test significance of fixed effects by backwards selection
modselect <- step(mod.all.full, lsmeans.calc=F, difflsmeans.calc=F, alpha.fixed=0.05)
modselect$anova.table
# Rebuild model omitting non-significant fixed effects
mod.all <- update(mod.all.full, formula(modselect$model))
# Identify outliers with standardized residuals > 2.5
out <- abs(residuals(mod.all)) > sd(residuals(mod.all)) * 2.5
Mcap.ff[out, ]  # outlying data points
# Refit model without outliers
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=Mcap.ff[!out, ])
# Print and save ANOVA table for model
anovatab <- anova(mod.all)
write.csv(round(anovatab, digits=3), file="anova_table.csv")
# pseudo-r2 value-- squared correlation between fitted and observed values
summary(lm(model.response(model.frame(mod.all)) ~ fitted(mod.all)))$r.squared

# Generate predictions and confidence intervals by parametric bootstrapping
pred.all <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")),
                         vis=factor(c("bleached", "not bleached")))
bootfit <- bootMer(mod.all, FUN=function(x) predict(x, pred.all, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
pred.all$fit <- predict(mod.all, pred.all, re.form=NA)
pred.all$lci <- apply(bootfit$t, 2, quantile, 0.05)
pred.all$uci <- apply(bootfit$t, 2, quantile, 0.95)

# Calculate when bleached are no longer different from healthy at each reef
# Test null hypotheses that S/H in bleached corals < non-bleached corals at each reef
#    on each day from 0 through 82, with p-value adjustments (method mvt) based on the 
#    number of tests on each day (3).
pvals <- matrix(0, nrow=83, ncol=3, dimnames=list(days=as.character(0:82), reef=levels(Mcap$reef)))
for(day in 0:82) {
  mod.ref.grid <- ref.grid(mod.all, at=list(days=day))
  pvals[day + 1, ] <- lsmeans::test(rbind(pairs(mod.ref.grid, by="reef")), adjust="mvt")$p.value
}
# Find which day/date pvalue is > 0.05 (no longer significantly different) at each reef
daystilrecov <- apply(pvals, 2, function(x) match(F, x < 0.05)) - 1
datestilrecov <- as.Date("2014-10-24", format="%F") + daystilrecov
datestilrecov; daystilrecov

#
# • Figure 5: Recovery dynamics --------------------------------------------------
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
      # Plot raw data +/- standard error
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
pdf("Figure5*.pdf", width=3.5, height=7)
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
# January (recovered)
Mcap.ff.jan <- subset(Mcap.ff, fdate=="20150114")
janmod <- lm(log(tot.SH) ~ reef * vis * tdom, data=Mcap.ff.jan)
anova(janmod)
janmean <- mean(log(Mcap.ff.jan$tot.SH))  # Grand geommean of january - no diff by reef or vis
exp(janmean)
# Comparisons
exp(octmeans[1])/exp(octmeans[2])  # Bleached are 88% lower than non-bleached in oct
exp(octmeans) / exp(janmean)  # Rel to jan, oct vals 95% lower in bleached, 61% lower in not-bleached
# May (post-recov)
Mcap.ff.may <- Mcap.ff[!out, ][which(Mcap.ff[!out, ]$fdate=="20150506"), ]
maymod <- lm(log(tot.SH) ~ reef * vis * tdom, data=Mcap.ff.may)
anova(maymod)
maymod <- lm(log(tot.SH) ~ vis, data=Mcap.ff.may)
anova(maymod)
data.frame(effect("vis", maymod))
exp(data.frame(effect("vis", maymod))$fit)
maymean <- mean(log(Mcap.ff.may$tot.SH))  # Grand geommean of january
exp(maymean)


# Nonbleached by date
mod <- lmer(log(tot.SH) ~ fdate + (1|colony), data=Mcap.ff[which(Mcap.ff$vis=="not bleached"), ])
anova(mod)
mod.ref.grid <- ref.grid(mod, at=list(fdate=c("20141024", "20150114", "20150506")))
contrast(mod.ref.grid, specs="fdate", method="pairwise")
lsmeans(mod.ref.grid, specs="fdate")

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
