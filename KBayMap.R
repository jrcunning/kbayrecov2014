library(RgoogleMaps)

# Import Kaneohe Bay map
KBay <- GetMap(center=c(21.460600, -157.809992), size=c(500,640),
               zoom=12, maptype="satellite", SCALE=2)
reef44 <- c(21.4767770, -157.8336070)
reef25 <- c(21.4611944, -157.8222500)
HIMB   <- c(21.4350000, -157.7910833)
reefcoords <- rbind(reef44, reef25, HIMB)

reefpts <- LatLon2XY.centered(KBay, reefcoords[,1], reefcoords[,2])


#pdf("KBay.pdf", width=7*(500/640), height=7)
PlotOnStaticMap(KBay)
points(reefpts$newX, reefpts$newY, col="yellow", cex=2, lwd=2)
text(reefpts$newX, reefpts$newY + 10, labels=c("Reef 44", "Reef 25", "HIMB"), 
     col="yellow", pos=3, cex=2)
#dev.off()


# Sampling dates stripchart
par(mar=c(2,4,2,4))
stripchart(c(0,11,31,53,82,194), pch=21, bg="black", bty="n", axes=F, cex=1.5)
axis(side=1, pos=1, xpd=T, labels=F, tcl=-0.5,
     at=as.numeric(as.Date(c("2014-10-01", "2014-11-01", "2014-12-01", "2015-01-01", "2015-02-01", 
                             "2015-03-01", "2015-04-01", "2015-05-01", "2015-06-01")) - as.Date("2014-10-24")))
axis(side=1, pos=1, xpd=T, tcl=0.3, lty=0,
     at=as.numeric(as.Date(c("2014-10-15", "2014-11-15", "2014-12-15", "2015-01-15", "2015-02-15", 
                             "2015-03-15", "2015-04-15", "2015-05-15")) - as.Date("2014-10-24")),
     labels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"))


# Outline map

install.packages("PBSmapping")
library(PBSmapping)
hawaii <- importGSHHS(file.path(.libPaths(), "PBSmapping/gshhs_f.b"), 
                              xlim=c(199, 208) , ylim=c(19, 24) , maxLevel=4)
plotMap(hawaii, col="grey", xlim=c(199,206), ylim=c(18,23))

oahu <- importGSHHS(file.path(.libPaths(), "PBSmapping/gshhs_f.b"), 
                     xlim=c(201.5, 202.5) , ylim=c(21, 22) , maxLevel=4)
plotMap(oahu, col="grey", xlim=c(201.7, 202.37), ylim=c(21.2, 21.8))


kbay <- importGSHHS(file.path(.libPaths(), "PBSmapping/gshhs_f.b"), 
                     xlim=c(202.1, 202.4) , ylim=c(21.4, 21.6) , maxLevel=4)
plotMap(kbay, col="grey", xlim=c(202.12, 202.28), ylim=c(21.4, 21.55))

install.packages("gmt")
library(gmt)
gmt(demo.par)
