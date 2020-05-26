##########################################
##### Figure_Generation.R ################
##### Generate Figures Used In Article ###
##### Kennedy-Shaffer and Lipsitch #######
##### Update: May 26, 2020 ###############
##### Please Cite https://doi.org/10.1101/2020.05.01.20087429
##### For questions, contact Lee Kennedy-Shaffer: lee_kennedyshaffer@g.harvard.edu
##########################################

### Packages Needed: ###
# install.packages("fields","dplyr")
require(fields)
require(dplyr)

# Parameter values needed from Outbreak_Simulations: direct_VE, R0


#### Begin from four combined data frames or matrices (one for PVals and one for VEs) of the results of:
## Analysis_IRTCRT_Res, Analysis_PH_Res, Analysis_MEMCPI_Res, Analysis_NPWP_Res, Analysis_SC_Res.
## Each row should be a simulation result with the columns as named in those outputs.
## Named VE.All.A and PVal.All.A (for SWT-A) and VE.All.B and PVal.All.B (for SWT-B). Only the SWT-A data frames need the IRT and CRT results.

Alpha <- 0.05

Power.All.A <- as.matrix(apply(X=PVal.All.A, MARGIN=2, FUN=function(x) mean(ifelse(x < Alpha,1,0), na.rm=TRUE)))
VE.Quantiles.A <- as.matrix(apply(X=VE.All.A, MARGIN=2, FUN=function(x) quantile(x, probs=c(0,.1,.25,.5,.75,.9,1),na.rm=TRUE)))
Power.All.B <- as.matrix(apply(X=PVal.All.B, MARGIN=2, FUN=function(x) mean(ifelse(x < Alpha,1,0), na.rm=TRUE)))
VE.Quantiles.B <- as.matrix(apply(X=VE.All.B, MARGIN=2, FUN=function(x) quantile(x, probs=c(0,.1,.25,.5,.75,.9,1),na.rm=TRUE)))


## Figure 1: Median and quartiles of estimates:
typenames <- c("iRCT","cRCT_gamma_coxph","PH","MEM.Ex1","CPI.Ex1",
               "lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
               "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
               "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")
typelabs <- c("IRT","CRT","SWT-PH","MEM","CPI","NPWP-1","NPWP-2","SC-1","SC-2","SC-Wt-1","SC-Wt-2")
ys <- seq(from=1,by=-.05,length.out=length(typenames))
if (direct_VE==0) {
  xlimval=c(-.6,.6)
} else {
  xlimval=c(0,1)
}

par(mfrow=c(1,2))
plot(x=VE.Quantiles.A["50%",typenames],
     y=ys, xlim=xlimval, ylim=c(min(ys)-.05,1),
     type="p", pch=18, cex=1.8, yaxt='n', ylab="",
     xlab="Estimated VE",
     main=bquote("A)"~"VE Estimates,"~R[0]==.(round(R0,2))*","~"VE"==.(direct_VE)*", SWT-A"))
for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
  abline(v=j, col="lightgray", lty="dotted")
}
for (i in 1:length(typenames)) {
  abline(h=ys[i], col="lightgray", lty="dotted")
}
abline(v=direct_VE, lty=5, col="gray", lwd=2)
points(x=VE.Quantiles.A["25%",typenames], y=ys, pch=20, cex=1.2)
points(x=VE.Quantiles.A["50%",typenames], y=ys, pch=18, cex=1.8)
points(x=VE.Quantiles.A["75%",typenames], y=ys, pch=20, cex=1.2)
for (i in 1:length(typenames)) {
  lines(x=VE.Quantiles.A[c("25%","75%"),typenames[i]], y=rep(ys[i],2), lty=1)
}
axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
       pt.cex=c(1.8,1.2), bg="white")

plot(x=VE.Quantiles.B["50%",typenames],
     y=ys, xlim=xlimval, ylim=c(min(ys)-.05,1),
     type="p", pch=18, cex=1.8, yaxt='n', ylab="",
     xlab="Estimated VE",
     main=bquote("B)"~"VE Estimates,"~R[0]==.(round(R0,2))*","~"VE"==.(direct_VE)*", SWT-B"))
for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
  abline(v=j, col="lightgray", lty="dotted")
}
for (i in 1:length(typenames)) {
  abline(h=ys[i], col="lightgray", lty="dotted")
}
abline(v=direct_VE, lty=5, col="gray", lwd=2)
points(x=VE.Quantiles.B["25%",typenames], y=ys, pch=20, cex=1.2)
points(x=VE.Quantiles.B["50%",typenames], y=ys, pch=18, cex=1.8)
points(x=VE.Quantiles.B["75%",typenames], y=ys, pch=20, cex=1.2)
for (i in 1:length(typenames)) {
  lines(x=VE.Quantiles.B[c("25%","75%"),typenames[i]], y=rep(ys[i],2), lty=1)
}
axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
       pt.cex=c(1.8,1.2), bg="white")
dev.off()


### Figure 2: Empirical Power or Type I Error:
xlimval=c(0,100)
typenames <- c("iRCT","cRCT_gamma_coxph","PH.Perm","MEM.Ex1","CPI.Ex1",
               "lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
               "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
               "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")
typelabs <- c("IRT","CRT","SWT-PH","MEM","CPI","NPWP-1","NPWP-2","SC-1","SC-2","SC-Wt-1","SC-Wt-2")
ys <- seq(from=1,by=-.05,length.out=length(typenames))
powerlab <- ifelse(direct_VE==0, "Type I Error", "Power")
plot(x=Power.All.A[typenames]*100,y=ys,
     xlim=xlimval, ylim=c(min(ys)-.05,1),
     type="p", pch=c(15,16,rep(24,8)), bg=rep(1,10), 
     yaxt='n', ylab="",
     xlab=paste0(powerlab, " (%)"),
     main=bquote(.(powerlab)~R[0]==.(round(R0,2))*","~"VE"==.(direct_VE)))
for (j in seq(from=xlimval[1],to=xlimval[2],by=10)) {
  abline(v=j, col="lightgray", lty="dotted")
}
for (i in 1:length(typenames)) {
  abline(h=ys[i], col="lightgray", lty="dotted")
}
points(x=Power.All.A[typenames]*100, y=ys,
       pch=c(15,16,rep(24,8)), bg=rep(1,10), cex=1.5)
points(x=Power.All.B[typenames[3:10]]*100, y=ys[3:10],
       pch=25, bg=1, cex=1.5)
axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
legend(x="bottom", legend=c("IRT","CRT","SWT-A","SWT-B"), pch=c(15,16,24,25), horiz=TRUE,
       pt.bg=rep(1,4), bg="white")


### Figures 3--5: Trends of Median Estimates and Power/Type I Error:
#### Begin with data frames or matrices named Medians.A, Medians.B, Power.A, and Power.B
##### Each data frame should have a column named R0 and the remaining columns named as in VE.All and Power.All above
##### Each row should have an R0 value and the median VE or power for the analysis types, under SWT-A or SWT-B
##### Only the A data frames need to have the IRT and CRT results.
if (direct_VE==0) {
  medianylim=c(-.5,.5)
  powerylim=c(0,100)
  medgrids=seq(-.4,.4,by=.2)
  powergrids=seq(0,100,by=20)
  ylab="Type I Error"
} else {
  medianylim=c(0,1)
  powerylim=c(0,100)
  medgrids=seq(0,1,by=.2)
  powergrids=seq(0,100,by=20)
  ylab="Power"
}

par(mfrow=c(1,2))
plot(x=Medians.A[,"R0"],y=Medians.A[,"iRCT"],ylim=medianylim,
     type="l",col=1,lty=1,
     xlab=bquote(R[0]),ylab=bquote("Median"~"VE Estimate"),
     main=bquote("A)"~"Median VE Estimate, VE"==.(direct_VE)))
points(x=Medians.A[,"R0"],y=Medians.A[,"iRCT"],pch=15)
lines(x=Medians.A[,"R0"],y=Medians.A[,"cRCT_gamma_coxph"],col=1,lty=2)
points(x=Medians.A[,"R0"],y=Medians.A[,"cRCT_gamma_coxph"],col=1,pch=16, cex=1.2)
lines(x=Medians.A[,"R0"],y=Medians.A[,"PH"],col=1,lty=3)
points(x=Medians.A[,"R0"],y=Medians.A[,"PH"],col=1,pch=23,bg=1)
lines(x=Medians.B[,"R0"],y=Medians.B[,"PH"],col=1,lty=3)
points(x=Medians.B[,"R0"],y=Medians.B[,"PH"],col=1,pch=23)
lines(x=Medians.A[,"R0"],y=Medians.A[,"lRR.NP.ToCaWt.Ex1"],col=1,lty=4)
points(x=Medians.A[,"R0"],y=Medians.A[,"lRR.NP.ToCaWt.Ex1"],col=1,pch=24,bg=1)
lines(x=Medians.B[,"R0"],y=Medians.B[,"lRR.NP.ToCaWt.Ex1"],col=1,lty=4)
points(x=Medians.B[,"R0"],y=Medians.B[,"lRR.NP.ToCaWt.Ex1"],col=1,pch=24)
lines(x=Medians.A[,"R0"],y=Medians.A[,"lRR.SC.ToCaWt.Ex1"],col=1,lty=5)
points(x=Medians.A[,"R0"],y=Medians.A[,"lRR.SC.ToCaWt.Ex1"],col=1,pch=25,bg=1)
lines(x=Medians.B[,"R0"],y=Medians.B[,"lRR.SC.ToCaWt.Ex1"],col=1,lty=5)
points(x=Medians.B[,"R0"],y=Medians.B[,"lRR.SC.ToCaWt.Ex1"],col=1,pch=25)
for (j in medgrids) {
  abline(h=j, col="lightgray", lty="dotted")
}
abline(h=VE, lty=5, col="gray", lwd=2)
legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A PH","SWT-A NPWP-2", "SWT-A SC-2","CRT","SWT-B PH","SWT-B NPWP-2","SWT-B SC-2"),
       lty=c(1,3,4,5,2,3,4,5),pch=c(15,23,24,25,16,23,24,25),pt.bg=c(NA,1,1,1,NA,NA,NA,NA), bg="white")

plot(x=Power.A[,"R0"],y=Power.A[,"iRCT"]*100,ylim=powerylim,type="l",col=1,lty=1,
     xlab=bquote(R[0]),ylab=bquote(.(ylab)~"(%)"),
     main=bquote("B)"~.(ylab)*", VE"==.(direct_VE)))
points(x=Power.A[,"R0"],y=Power.A[,"iRCT"]*100,pch=15)
lines(x=Power.A[,"R0"],y=Power.A[,"cRCT_gamma_coxph"]*100,col=1,lty=2)
points(x=Power.A[,"R0"],y=Power.A[,"cRCT_gamma_coxph"]*100,col=1,pch=16, cex=1.2)
lines(x=Power.A[,"R0"],y=Power.A[,"PH.Perm"]*100,col=1,lty=3)
points(x=Power.A[,"R0"],y=Power.A[,"PH.Perm"]*100,col=1,pch=23,bg=1)
lines(x=Power.B[,"R0"],y=Power.B[,"PH.Perm"]*100,col=1,lty=3)
points(x=Power.B[,"R0"],y=Power.B[,"PH.Perm"]*100,col=1,pch=23)
lines(x=Power.A[,"R0"],y=Power.A[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,lty=4)
points(x=Power.A[,"R0"],y=Power.A[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,pch=24,bg=1)
lines(x=Power.B[,"R0"],y=Power.B[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,lty=4)
points(x=Power.B[,"R0"],y=Power.B[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,pch=24)
lines(x=Power.A[,"R0"],y=Power.A[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,lty=5)
points(x=Power.A[,"R0"],y=Power.A[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,pch=25,bg=1)
lines(x=Power.B[,"R0"],y=Power.B[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,lty=5)
points(x=Power.B[,"R0"],y=Power.B[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,pch=25)
for (j in powergrids) {
  abline(h=j, col="lightgray", lty="dotted")
}
if (VE==0) {
  abline(h=5, lty=5, col="gray", lwd=2)
}
dev.off()


### Figure 6: Average Period-Specific Estimates:
#### Begin with an array where the first two dimensions are the row-binded results of 
##### Analysis_NPWP_ByPeriod and Analysis_SC_ByPeriod, named Analysis_ByPeriod.A (for SWT-A) and Analysis_ByPeriod.B (for SWT-B)
hm.colors <- snow.colors(n=72)[65:2]
hm.colors.gray <- gray.colors(n=64, start=.8, end=.2)

xlimval.A <- c(2.5,10.5)
xseq.A <- seq(.5,11.5,by=1)
xlabels.A <- 3:10
xlimval.B <- c(5.5,42.5)
xseq.B <- seq(0.5,44.5,by=1)
xlabels.B <- seq(6,42,by=2)

Heat.A <- 1-exp(apply(Analysis_ByPeriod.A[c("SCWt.lRR.Ex1","SC.lRR.Ex1","NP.lRR.Ex1"),,],
                     c(1,2), mean, na.rm=TRUE))
Heat.A <- ifelse(Heat.A < -.2, -.2, Heat.A) # Restricts to -.2 so that differences can be perceived
Heat.B <- 1-exp(apply(Analysis_ByPeriod.B[c("SCWt.lRR.Ex1","SC.lRR.Ex1","NP.lRR.Ex1"),,],
                     c(1,2), mean, na.rm=TRUE))
Heat.B <- ifelse(Heat.B < -.2, -.2, Heat.B)
par(mfrow=c(2,1))
image.plot(x=xseq.A,y=1:4,z=t(Heat.A),xlim=xlimval.A,ylim=c(1,4),
           ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
           main=bquote("A)"~"Estimated VE by Period, SWT-A"),
           col=hm.colors)
axis(side=2,
     labels=c("SC-Wt","SC","NPWP"),
     at=c(1.5,2.5,3.5))
axis(side=1, labels=xlabels.A, at=xlabels.A)
abline(h=2, col=1, lty=1, lwd=1.5)
abline(h=3, col=1, lty=1, lwd=1.5)

image.plot(x=xseq.B,y=1:4,z=t(Heat.B),xlim=xlimval.B,ylim=c(1,4),
           ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
           main=bquote("B)"~"Estimated VE by Period, SWT-B"),
           col=hm.colors)
axis(side=2,
     labels=c("SC-Wt","SC","NPWP"),
     at=c(1.5,2.5,3.5))
axis(side=1, labels=xlabels.B, at=xlabels.B)
abline(h=2, col=1, lty=1, lwd=1.5)
abline(h=3, col=1, lty=1, lwd=1.5)
dev.off()
