source("~/Documents/Academia/HIMB/steponeR/steponeR.R")
files <- list("20150910_KBayRecov_Mcap_IntraColVar_data.csv", 
              "20150910_KBayRecov_Mcap_IntraColVar_2_data.csv")

icv <- steponeR(files, target.ratios=c("C.Mcap", "D.Mcap"),
                fluor.norm=list(C=2.26827, D=0, Mcap=0.84815),
                copy.number=list(C=33, D=3, Mcap=1),
                ploidy=list(C=1, D=1, Mcap=2),
                extract=list(C=0.813, D=0.813, Mcap=0.982))
icv$result$C.Mcap[which(is.na(icv$result$C.Mcap))] <- 0
icv$result$D.Mcap[which(is.na(icv$result$D.Mcap))] <- 0
icv$result$tot.SH <- icv$result$C.Mcap + icv$result$D.Mcap

icv$result$colony <- cut(as.numeric(icv$result$Sample.Name), include.lowest=TRUE,
                  breaks=c(701, 706, 712, 718, 724, 730, 736), 
                  labels=c("31", "32", "25", "44", "20", "3"))

boxplot(tot.SH ~ colony, data=icv$result)
boxplot(log(tot.SH) ~ colony, data=icv$result, xlab="Colony", ylab="ln S/H")
mod <- aov(log(tot.SH) ~ colony, data=icv$result)
summary(mod)
