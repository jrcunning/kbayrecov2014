# Title: Dynamics and recovery of Montipora capitata symbioses following bleaching in Kaneohe Bay
# Author: Ross Cunning
# Last updated: 27 August, 2015

source("setup.R")
source("gsp.R")

# =================================================================================================
# DATA ANALYSIS
# =================================================================================================
# • Patterns in dominant symbionts and bleaching --------------------------------------------
# --Use Mcap.f data.frame (duplicates removed but data unsuitable for quantitative analysis remain)
# --Figure: Symbiont clades in each colony over time -------------------------
clades <- melt(Mcap.f, id.vars=c("sample", "date", "vis", "reef"), measure.vars="syms",
               factorsAsStrings=FALSE)
plot(value ~ date, clades[clades$vis=="bleached", ])
plot(value ~ date, clades[clades$vis=="not bleached", ])
# Create matrix for image function
clades$value <- as.numeric(factor(clades$value))
clades <- dcast(clades, vis + sample + reef ~ date, drop=T)
clades[is.na(clades)] <- -1  # Recode missing values as -1
clades <- clades[with(clades, order(rev(vis), clades[, 4], clades[, 5], clades[, 6], 
                                    clades[, 7], clades[, 8], clades[, 9])), ]
clades.m <- as.matrix(clades[,4:9], row.names=as.character(clades$sample))
rownames(clades.m) <- as.character(clades$sample)
# Plot heatmap image
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
# --Figure: D:C ratio in each colony over time --------------------------
logDC <- melt(Mcap.f, id.vars=c("sample", "date", "vis", "reef"), measure.vars=c("logDC"),
                 factorsAsStrings=FALSE)
logDC$value[which(logDC$valu==-Inf)] <- -12
logDC$value[which(logDC$valu==Inf)] <- 12
head(logDC)
xyplot(value ~ date | reef, groups = ~ sample, data = logDC, type="o")
logDC <- dcast(logDC, reef + vis + sample ~ date, drop=T)
logDC <- logDC[with(logDC, order(vis, reef, sample)), ]

rownames(logDC) <- logDC$sample
logDC <- logDC[with(logDC, order(rev(vis), logDC[,4], logDC[,5], logDC[,6], logDC[,7],
                                 logDC[,8], logDC[,9])), ]
logDC.m <- as.matrix(logDC[,4:9])
rownames(logDC.m) <- rownames(logDC)
# Find range of non-infinite ratios
range(logDC.m[which(!is.infinite(logDC.m))], na.rm=T)
# Assign -Inf=-12, +Inf=+12
logDC.m[which(logDC.m==-Inf)] <- -12
logDC.m[which(logDC.m==Inf)] <- 12

# Plot D:C ratio heatmap
par(mar=c(5,7,3,1), bg="white")
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
# --Statistical tests of symbiont clade, reef, and visual status--------------------
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
# • Quantitative analysis of S/H ratios-- data setup ------------------------------------------------------------
# Remove CD colony 40?


# Bleaching severity in October -- effect of symbiont clade and visual appearance
Mcap.ff.oct <- Mcap.ff[which(Mcap.ff$fdate=="20141024"), ]
model <- lm(log(tot.SH) ~ dom:vis, data=Mcap.ff.oct)
anova(model, test="F")  # tdom:vis significant, (reef is not significant if included)
plot(effect("dom:vis", model))
TukeyHSD(aov(model))  # BLEACHED C SAME AS NONBLEACHED D--- BECAUSE OF COPY NUMBER?
# Conclusion: Bleached C colonies had lowest S/H, Unbleached C and D colonies had same S/H



### MIXED MODELS TO ANALYZE RECOVERY TRENDS
# Subset data for mixed models
Mcap.ff.octdec <- subset(Mcap.ff, date <= "2014-12-16")
Mcap.ff.bleached <- subset(Mcap.ff, vis=="bleached")
Mcap.ff.bleached.octdec <- subset(Mcap.ff.bleached, date <= "2014-12-16")
Mcap.ff.bleached.octjan <- subset(Mcap.ff.bleached, date <= "2015-01-14")
Mcap.ff.notbleached <- subset(Mcap.ff, vis=="not bleached")
Mcap.ff.notbleached.octdec <- subset(Mcap.ff.notbleached, date <= "2014-12-16")
Mcap.ff.notbleached.octjan <- subset(Mcap.ff.notbleached, date <= "2015-01-14")

# HOW DID BLEACHED CORALS RECOVER OVER TIME AT DIFFERENT REEFS?
# Visualize data/ trajectories for each colony over time
xyplot(log(tot.SH) ~ days | reef, groups = ~ sample, data = Mcap.ff.bleached[order(Mcap.ff.bleached$days), ],
       type = "o", layout=c(3, 1), main="bleached colonies")

# ---first, treat date as factor, not continuous ------------------
model <- lmer(log(tot.SH) ~ fdate + reef + fdate:reef + (1|sample), data=Mcap.ff.bleached)
dropterm(model, test="Chisq")   # significant date:reef interaction with all time points
                                # smaller p value with unfiltered data
plot(effect("fdate:reef", model), multiline=T, ci.style="bars")

model.lsm <- lsmeans(model, c("fdate", "reef"), contr = "cld")
model.lsm   # pairwise tests of means

model.lsm$date <- as.Date(model.lsm$fdate, format="%Y%m%d")
model.lsm <- model.lsm[order(model.lsm$date), ]
par(mfrow=c(1,1))
plot(lsmean ~ date, model.lsm[which(model.lsm$reef=="HIMB"), ], col="red", type="o", ylim=c(-3,2),
     xlab="Date", ylab="log S/H ratio", main="Recovery of bleached corals")
lines(lsmean ~ date, model.lsm[which(model.lsm$reef=="25"), ], col="darkgreen", type="o")
lines(lsmean ~ date, model.lsm[which(model.lsm$reef=="44"), ], col="blue", type="o")
arrows(model.lsm[which(model.lsm$reef=="44"), ]$date, 
       model.lsm[which(model.lsm$reef=="44"), ]$lsmean+model.lsm[which(model.lsm$reef=="44"), ]$SE, 
       model.lsm[which(model.lsm$reef=="44"), ]$date, 
       model.lsm[which(model.lsm$reef=="44"), ]$lsmean-model.lsm[which(model.lsm$reef=="44"), ]$SE,
       code=0, col="blue")
arrows(model.lsm[which(model.lsm$reef=="25"), ]$date, 
       model.lsm[which(model.lsm$reef=="25"), ]$lsmean+model.lsm[which(model.lsm$reef=="25"), ]$SE, 
       model.lsm[which(model.lsm$reef=="25"), ]$date, 
       model.lsm[which(model.lsm$reef=="25"), ]$lsmean-model.lsm[which(model.lsm$reef=="25"), ]$SE,
       code=0, col="darkgreen")
arrows(model.lsm[which(model.lsm$reef=="HIMB"), ]$date, 
       model.lsm[which(model.lsm$reef=="HIMB"), ]$lsmean+model.lsm[which(model.lsm$reef=="HIMB"), ]$SE, 
       model.lsm[which(model.lsm$reef=="HIMB"), ]$date, 
       model.lsm[which(model.lsm$reef=="HIMB"), ]$lsmean-model.lsm[which(model.lsm$reef=="HIMB"), ]$SE,
       code=0, col="red")
legend("topright", levels(model.lsm$reef)[c(2,1,3)], lty=1, col=c("blue", "darkgreen", "red"))

xyplot(lsmean ~ date, groups = ~reef, data=model.lsm[order(model.lsm$date), ], type="o",
       auto.key=T)

# ---treat time as continuous, model linear response from october until december------------------
xyplot(log(tot.SH) ~ days | reef, groups = ~ sample, data = Mcap.ff.bleached.octdec,
       type = "o", layout=c(3, 1), main="bleached colonies", auto.key=F)
model <- lmer(log(tot.SH) ~ days * reef + (days|sample), data=Mcap.ff.bleached.octdec)  # random slopes and intercepts for each colony (sample)
sel <- lmerTest::step(model); sel  # Backward elimination of non-significant random and fixed effects
model <- update(model, formula=formula(sel$model))  # Update model to include only significant terms
dropterm(model, test="Chisq")  # LRT for fixed effects
# Model diagnostic plots
plot(model); qqmath(model)
1-var(residuals(model))/(var(model.response(model.frame(model))))  #pseudo-R^2
# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(model, Mcap.ff.bleached.octdec, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]
model <- update(model, data = rm.outliers$data)
plot(effect("days:reef", model, KR=T), multiline=T, ci.style="bands")
model.lst <- lstrends(model, ~ reef, var = "days")
cld(model.lst)  # SLOPES ARE SIG. DIFF BETWEEN HIMB AND 44, 25 IS N.S. INTERMEDIATE

# ---POLYNOMIAL REGRESSION OVER RECOVERY PERIOD-----------
# Plot raw data for each colony
xyplot(log(tot.SH) ~ days | reef, groups = ~ sample, data = Mcap.ff.bleached.octjan,
       type = "o", layout=c(3, 1), main="bleached colonies", auto.key=F)
# Build model
lmodel <- lme4::lmer(log(tot.SH) ~ days * reef + (days|sample), data=Mcap.ff.bleached.octjan)
qmodel <- lme4::lmer(log(tot.SH) ~ poly(days, 2) * reef + (poly(days,2)|sample), data=Mcap.ff.bleached.octjan)
cmodel <- lme4::lmer(log(tot.SH) ~ poly(days, 3) * reef + (poly(days,3)|sample), data=Mcap.ff.bleached.octjan)
anova(lmodel, qmodel, cmodel)  # quadratic model is best
qmodel2 <- lme4::lmer(log(tot.SH) ~ poly(days, 2) * reef + (1|sample), data=Mcap.ff.bleached.octjan)
anova(qmodel, qmodel2)  # model with random intercept only is best
model <- qmodel2
dropterm(model, test="Chisq")  # LRT indicates interaction term is significant (compared to model without interaction)
summary(model)  # Model parameters
mcp.fnc(model)  # Model diagnostics -- residuals are homoscedastic, QQplot is normal

# # Remove outliers with residuals > 2.5 s.d.'s from 0
# rm.outliers <- romr.fnc(model, Mcap.ff.bleached.octjan, trim=2.5)
# rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
# model <- update(model, data = rm.outliers$data)

# Make predictions
# Predictions using base R
pred <- expand.grid(days=seq(0,82,2), reef=factor(c("44", "25", "HIMB")))
pred$fit <- predict(qmodel, newdata=pred, re.form=NA, type="response")  # re.form=NA means all random effect values are set to zero--predictions are made at the population level
xyplot(fit ~ days, groups=reef, data=pred, type="l")  # Plot model predictions

# Generate predictions and prediction intervals using bootMer over days 0-82
newdat <- expand.grid(days=seq(0,82,2), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(model, FUN=function(x) predict(x, newdat, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat$pred <- predict(model, newdat, re.form=NA)
newdat$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat$uci <- apply(bootfit$t, 2, quantile, 0.95)
# Plot predicted values with 90% confidence intervals
with(newdat[newdat$reef=="HIMB", ], {
  plot(pred ~ days, type="l", col="red", ylim=c(-3.25,1.5))
  addpoly(days, lci, uci, col=alpha("red", 0.5))
})
with(newdat[newdat$reef=="25", ], {
  lines(pred ~ days, col="blue")
  addpoly(days, lci, uci, col=alpha("blue", 0.5))
})
with(newdat[newdat$reef=="44", ], {
  lines(pred ~ days, col="darkgreen")
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.5))
})


# Generate prediction intervals using effects package
eff <- Effect(c("days", "reef"), model, xlevels=list(days=1:82), confidence.level=0.90)
#plot(eff, multiline=T, ci.style="bands", ylim=c(-3.5,2))
effdf <- data.frame(eff)   # HOW TO INTERPRET SE VS. CONF.INT.?
with(effdf[effdf$reef=="HIMB", ], {
  plot(fit ~ days, type="l", col="red", ylim=c(-3.5,2))
  addpoly(days, lower, upper, col=alpha("red", 0.5))
})
with(effdf[effdf$reef=="25", ], {
  lines(fit ~ days, col="blue")
  addpoly(days, lower, upper, col=alpha("blue", 0.5))
})
with(effdf[effdf$reef=="44", ], {
  lines(fit ~ days, col="darkgreen")
  addpoly(days, lower, upper, col=alpha("darkgreen", 0.5))
})

# ---PIECEWISE POLYNOMIAL REGRESSION OVER FULL TIME SERIES-----------------
# Fit model--------------------
# Spline fit with knot at 82 days (january), quadratic fit before then, linear after. continuous at knot.
sp <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
model <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (1 | sample), data=Mcap.ff.bleached)
model2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (sp(days) | sample), data=Mcap.ff.bleached)
anova(model, model2)  # Model with random intercept only is best
summary(model)

# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(model, Mcap.ff.bleached, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
model <- update(model, data = rm.outliers$data)

# Generate predictions and prediction intervals using bootMer over days 0-194
newdat <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(model, FUN=function(x) predict(x, newdat, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat$pred <- predict(model, newdat, re.form=NA)
newdat$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat$uci <- apply(bootfit$t, 2, quantile, 0.95)
newdat <- split(newdat, f=newdat$reef)
# Summary data
mdf <- model.frame(model)
model.data.summ <- data.frame(expand.grid(reef=levels(mdf$reef), days=as.numeric(as.character(levels(factor(mdf$`sp(days)`[,1]))))),
                              mean=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp(days)`[,1])), FUN=sd)$x)

model.data.summ <- split(model.data.summ, f=model.data.summ$reef)
# Plot raw data mean ± sd with model fit and CI-------------

par(mfrow=c(2,2), mar=c(2,2,2,2))
with(model.data.summ$`HIMB`, {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.25,2.5), bty="n")
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat$"HIMB", lines(days, pred))
  with(newdat$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
with(model.data.summ$`25`, {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.25,2.5), bty="n")
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat$"25", lines(days, pred))
  with(newdat$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
with(model.data.summ$`44`, {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.25,2.5), bty="n")
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat$"44", lines(days, pred))
  with(newdat$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(2,1,1,2))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
box(lty=2, col=alpha("black", 0.8))
with(newdat$"HIMB", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat$"25", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat$"44", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# Plot same data with different layout-------------
#layout(mat = rbind(c(1,4),c(2,4),c(3,5)))
layout(mat=matrix(c(1,2,3,4,4)))
#par(mfrow=c(3,1))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1))
par(mar=c(0,2,0,0))
with(model.data.summ$`44`, {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat$"44", lines(days, pred))
  with(newdat$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`25`, {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat$"25", lines(days, pred))
  with(newdat$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`HIMB`, {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat$"HIMB", lines(days, pred))
  with(newdat$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat$"HIMB", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat$"25", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat$"44", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# # Add temperature plot
# par(mar=c(1,1,0,2))
# with(rf44[which(rf44$date >= "2014-10-24" & rf44$date <= "2015-01-14"), ], {
#   plot(mean ~ date, type="l", col="darkgreen", ylim=c(22,29), xaxs="i", yaxs="i")
#   addpoly(date, min, max, col=alpha("darkgreen", 0.3))
# })
# with(rf25[which(rf25$date >= "2014-10-24" & rf25$date <= "2015-01-14"), ], {
#   points(mean ~ date, type="l", col="blue", ylim=c(22,29))
#   addpoly(date, min, max, col=alpha("blue", 0.3))
# })
# with(HIMB[which(HIMB$date >= "2014-10-24" & HIMB$date <= "2015-01-14"), ], {
#   points(mean ~ date, type="l", col="red", ylim=c(22,29))
#   addpoly(date, min, max, col=alpha("red", 0.3))
# })

















# Look at not-bleached colonies----------------
# remove CD colony
Mcap.ff.notbleached <- droplevels(subset(Mcap.ff.notbleached, tdom!="CD"))
Mcap.ff.notbleached.octdec <- subset(Mcap.ff.notbleached, date <= "2014-12-16")
Mcap.ff.notbleached.octjan <- subset(Mcap.ff.notbleached, date <= "2015-01-14")
xyplot(log(tot.SH) ~ days | reef + tdom, groups = ~ sample, data=Mcap.ff.notbleached[order(Mcap.ff.notbleached$days), ],
       type = "o", layout=c(3, 2), main="not bleached colonies")
# piecewise polynomial regression
sp <- function(x) gsp(x, knots=c(82), degree=c(2,1), smooth=0)
model <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (1 | sample), data=Mcap.ff.notbleached)
model2 <- lmerTest::lmer(log(tot.SH) ~ sp(days) * reef + (sp(days) | sample), data=Mcap.ff.notbleached)
anova(model, model2)  # Model with random intercept only is best
sp2 <- function(x) gsp(x, knots=c(82), degree=c(1,1), smooth=0)
model3 <- lmerTest::lmer(log(tot.SH) ~ sp2(days) * reef + (1 | sample), data=Mcap.ff.notbleached)
anova(model, model3)
model <- model3
summary(model)
# Remove outliers with residuals > 2.5 s.d.'s from 0
rm.outliers <- romr.fnc(model, Mcap.ff.notbleached, trim=2.5)
rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2.5), ]  # which data points are removed
model <- update(model, data = rm.outliers$data)
# Generate predictions and prediction intervals using bootMer over days 0-194
newdat <- expand.grid(days=seq(0,194,1), reef=factor(c("44", "25", "HIMB")))
bootfit <- bootMer(model, FUN=function(x) predict(x, newdat, re.form=NA), nsim=999)
# Extract 90% confidence interval on predicted values
newdat$pred <- predict(model, newdat, re.form=NA)
newdat$lci <- apply(bootfit$t, 2, quantile, 0.05)
newdat$uci <- apply(bootfit$t, 2, quantile, 0.95)
newdat <- split(newdat, f=newdat$reef)
# Plot same data with different layout-------------
mdf <- model.frame(model)
model.data.summ <- data.frame(expand.grid(reef=levels(mdf$reef), days=as.numeric(as.character(levels(factor(mdf$`sp2(days)`[,1]))))),
                              mean=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp2(days)`[,1])), FUN=mean)$x,
                              sd=aggregate(mdf$`log(tot.SH)`, by=list(interaction(mdf$reef, mdf$`sp2(days)`[,1])), FUN=sd)$x)

model.data.summ <- split(model.data.summ, f=model.data.summ$reef)
#layout(mat = rbind(c(1,4),c(2,4),c(3,5)))
layout(mat=matrix(c(1,2,3,4,4)))
#par(mfrow=c(3,1))
par(mgp=c(2,0.4,0), oma=c(1,1,1,1))
par(mar=c(0,2,0,0))
with(model.data.summ$`44`, {
  plot(mean ~ days, pch=21, bg="darkgreen", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05)
  with(newdat$"44", lines(days, pred))
  with(newdat$"44", addpoly(days, lci, uci, col=alpha("darkgreen", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`25`, {
  plot(mean ~ days, pch=21, bg="blue", ylim=c(-4.1,2.4), bty="n", xaxt="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat$"25", lines(days, pred))
  with(newdat$"25", addpoly(days, lci, uci, col=alpha("blue", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(0,2,0,0))
with(model.data.summ$`HIMB`, {
  plot(mean ~ days, pch=21, bg="red", ylim=c(-4.1,2.4), bty="n", tck=-0.03)
  arrows(days, mean+sd, days, mean-sd, code=3, angle=90, length=0.05, xpd=T)
  with(newdat$"HIMB", lines(days, pred))
  with(newdat$"HIMB", addpoly(days, lci, uci, col=alpha("red", 0.3)))
  rect(xleft = 0, ybottom = -3, xright = 82, ytop = 1, lty = 2, border=alpha("black", 0.8))
})
par(mar=c(1,2,5,0))
plot(NA, ylim=c(-3,1), xlim=c(0,82), xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
mtext(side=3, text = "Days", line=2.5, cex=0.75)
box(lty=2, col=alpha("black", 0.8))
with(newdat$"HIMB", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("red", 0.3))
})
with(newdat$"25", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("blue", 0.3))
})
with(newdat$"44", {
  lines(days, pred)
  addpoly(days, lci, uci, col=alpha("darkgreen", 0.3))
})
# # Add temperature plot
# par(mar=c(1,1,0,2))
# with(rf44[which(rf44$date >= "2014-10-24" & rf44$date <= "2015-01-14"), ], {
#   plot(mean ~ date, type="l", col="darkgreen", ylim=c(22,29), xaxs="i", yaxs="i")
#   addpoly(date, min, max, col=alpha("darkgreen", 0.3))
# })
# with(rf25[which(rf25$date >= "2014-10-24" & rf25$date <= "2015-01-14"), ], {
#   points(mean ~ date, type="l", col="blue", ylim=c(22,29))
#   addpoly(date, min, max, col=alpha("blue", 0.3))
# })
# with(HIMB[which(HIMB$date >= "2014-10-24" & HIMB$date <= "2015-01-14"), ], {
#   points(mean ~ date, type="l", col="red", ylim=c(22,29))
#   addpoly(date, min, max, col=alpha("red", 0.3))
# })







# model not bleached colonies oct - dec------------
model <- lmer(log(tot.SH) ~ days * tdom + reef + (1+days|sample), data=Mcap.ff.notbleached.octdec)
dropterm(model, test="Chisq")
plot(effect("days:tdom", model))
model.lst <- lstrends(model, "tdom", var="days")
cld(model.lst)  # TRENDS ARE SIG DIFF, C GOES DOWN< D GOES UP







# model not bleached colonies oct - jan
model <- lmer(log(tot.SH) ~ poly(days, 2) * tdom + reef + (1+days|sample), data=Mcap.ff.notbleached.octjan)
#mcp.fnc(model)
#rm.outliers <- romr.fnc(model, Mcap.ff.notbleached.octjan)
#rm.outliers$data0[which(abs(rm.outliers$data0$rstand) > 2), ]
#model.data <- rm.outliers$data
#model <- update(model, data = model.data)


dropterm(model, test="Chisq")
plot(Effect(c("days", "tdom"), model), multiline=T, ci.style="bands")








# Look at presence of C and D in bleached corals at each time point
Mcap.f.bleached$present <- ifelse(is.na(Mcap.f.bleached$C.SH), ifelse(is.na(Mcap.f.bleached$D.SH), NA, "D"),
                                  ifelse(is.na(Mcap.f.bleached$D.SH), "C", "CD"))
plot(factor(present) ~ factor(date), Mcap.f.bleached)
plot(factor(present) ~ factor(reef), Mcap.f.bleached)
plot(factor(present) ~ interaction(factor(date), factor(reef)), Mcap.f.bleached)


# Look at abundance of D in bleached corals
boxplot(log(Mcap.f.bleached$D.SH) ~ Mcap.f.bleached$days)
boxplot(log(Mcap.f.bleached$D.SH) ~ Mcap.f.bleached$reef)

model <- lmer(log(D.SH) ~ fdate * reef + (1|sample), data=Mcap.f.bleached)
summary(model)
dropterm(model, test="Chisq")
lsmeans(model, c("fdate", "reef"), contr="cld")
boxplot(log(D.SH) ~ fdate * reef, Mcap.f.bleached)


xyplot(log(D.SH) ~ days | reef, group = ~sample, data=Mcap.f.bleached[order(Mcap.f.bleached$days), ], type="o")
xyplot(log(D.SH) ~ days | reef + tdom, group = ~sample, data=Mcap.f.notbleached[order(Mcap.f.notbleached$days), ], type="o")


xyplot(log(D.SH) ~ days | sample, data=Mcap.f.C[order(Mcap.f.C$days), ], type="o")

xyplot(log(C.SH) ~ days | sample, data=Mcap.f.D[order(Mcap.f.D$days), ], type="o")

Mcap.f.bleached$DC <- Mcap.f.bleached$D.SH / Mcap.f.bleached$C.SH


xyplot(DC ~ days | sample, data=Mcap.f.bleached[order(Mcap.f.bleached$days), ], type="o")





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
