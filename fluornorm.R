

file <- "20150907_KBayRecov_Mcap_fluornorm_data.csv"
fn <- read.csv(file, skip=6, na.strings="Undetermined")[, c(1:3, 7, 10)]
fn <- droplevels(fn[which(fn$Sample.Name!=""), ])

# Plot standard curves
plot(C_ ~ log10(Quantity), data=fn, pch=21, bg=fn$Target.Name, ylab="CT value")
legend("topright", levels(fn$Target.Name), pch=21, pt.bg=factor(levels(fn$Target.Name)))

# Sets clade D target as baseline contrast
contrasts(fn$Target.Name) <- contr.treatment(levels(fn$Target.Name), base=2)

mod <- lm(C_ ~ log10(Quantity) + Target.Name, data=fn)
anova(mod)
summary(mod)
coef(mod)  # Fluorescence normalization values: 0 for clade D, 0.84815 for Mcap, and 2.26827 for clade C

# Plot fitted values for each target
newdat <- expand.grid(Quantity=10^seq(2,6,1), Target.Name=factor(c("D", "C", "Mcap")))
newdat$fit <- predict(mod, newdat)
with(newdat[which(newdat$Target.Name=="D"), ], lines(log10(Quantity), fit, col="red"))
with(newdat[which(newdat$Target.Name=="Mcap"), ], lines(log10(Quantity), fit, col="green"))
with(newdat[which(newdat$Target.Name=="C"), ], lines(log10(Quantity), fit, col="black"))
