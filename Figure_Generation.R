
basefolder <- "/Users/leekennedy-shaffer/Documents/Harvard/Research/ID Time-Varying Estimand/Code/Outbreak_Sims"
Compfolder <- paste0(basefolder,"/Compiled")
figfolder <- paste0(basefolder,"/Figures")
load(paste0(basefolder,"/Varying.Rda"))

NumSims <- 500
scens <- 3:4
params <- c("01","02","03","04","05","06","07","08","09","10","11","12")

### Packages Needed: ###
library(fields)
library(dplyr)

VE.Medians <- NULL
PVal.All <- NULL
Power.All <- NULL

for (scen in scens) {
  for (param in params) {
    name <- paste0("S",scen,"_P",param)
    paramRow <- Varying[((scen-1)*12+as.numeric(param)),]
    
    ### Load all Result Data Sets: ###
    load(paste0(Compfolder,"/Compiled_",name,".Rda"))
    load(paste0(Compfolder,"/MEMEsts_",name,".Rda"))
    load(paste0(Compfolder,"/MEMPVals_",name,".Rda"))
    load(paste0(Compfolder,"/TERD_",name,".Rda"))
    load(paste0(Compfolder,"/TERR_",name,".Rda"))
    load(paste0(Compfolder,"/TElRR_",name,".Rda"))
    load(paste0(Compfolder,"/PVals_NP_",name,".Rda"))
    load(paste0(Compfolder,"/FullRes2Comp_",name,".Rda"))
    if (scen==3 | param %in% c("02","03","05","06","07","08","09","11","12")) {
      load(paste0(Compfolder,"/PVals_SC_",name,".Rda"))
    }
    
    VE.SP <- cbind(data.frame(Sim=1:NumSims), 
                 get(paste0("Compiled_",name))$VEs[,c("iRCT","cRCT_gaussian_coxme","cRCT_gaussian_coxph","cRCT_gamma_coxph","cRCT_gee")])
    PVal.SP <- cbind(data.frame(Sim=1:NumSims), 
                   get(paste0("Compiled_",name))$pvals[,c("iRCT","cRCT_gaussian_coxme","cRCT_gaussian_coxph","cRCT_gamma_coxph","cRCT_gee")])
    VEs.MEM <- Ests
    VEs.MEM[,-1] <- 1-exp(VEs.MEM[,-1])
    VE.SP <- merge(VE.SP, VEs.MEM, by="Sim")
    PVal.SP <- merge(PVal.SP, PVals, by="Sim")
    for (type in c("RD","RR","lRR")) {
      base <- get(paste0("TE.",type))
      names <- colnames(base)
      colnames(base) <- c("Sim",paste0(type,".",names[-1]))
      if (type=="RD") {
        base[,-1] <- -1*base[,-1]
      } else if (type=="RR") {
        base[,-1] <- 1-base[,-1]
      } else if (type=="lRR") {
        base[,-1] <- 1-exp(base[,-1])
      }
      VE.SP <- merge(VE.SP, base, by="Sim")
    }
    PVal.SP <- merge(PVal.SP, NPPVals, by="Sim")
    if (scen==3 | param %in% c("02","03","05","06","07","08","09","11","12")) {
      PVal.SP <- merge(PVal.SP, SC_PVals, by="Sim")
    }
    Power.SP <- apply(PVal.SP[,-1], 2, function(x) mean(ifelse(x < 0.05,1,0), na.rm=TRUE))
    
    Res.SP <- apply(VE.SP, 2, mean, na.rm=TRUE)
    Res.SP <- rbind(Res.SP, apply(VE.SP, 2, sd, na.rm=TRUE))
    rownames(Res.SP) <- c("Mean","SD")
    Res.SP <- rbind(Res.SP, apply(VE.SP, 2, 
                            function(x) quantile(x, probs=c(0,.1,.25,.5,.75,.9,1), 
                                                 na.rm=TRUE)))
    
    VE.Medians <- bind_rows(VE.Medians, c(Scen=scen,Param=as.numeric(param),
                                      Beta=paramRow$Beta,R0=paramRow$R0,VE=paramRow$VE,
                                      Res.SP["50%",-1]))
    PVal.All <- bind_rows(PVal.All, cbind(data.frame(Scen=rep(scen,dim(PVal.SP)[1]),Param=rep(as.numeric(param),dim(PVal.SP)[1]),
                                                 Beta=rep(paramRow$Beta,dim(PVal.SP)[1]),R0=rep(paramRow$R0,dim(PVal.SP)[1]),
                                                 VE=rep(paramRow$VE,dim(PVal.SP)[1])),
                                      PVal.SP))
    Power.All <- bind_rows(Power.All, c(Scen=scen,Param=as.numeric(param),
                                    Beta=paramRow$Beta,R0=paramRow$R0,VE=paramRow$VE,
                                    Power.SP))
    assign(x=paste0("fullResComp_",name), value=fullResComp)
    assign(x=paste0("VE_",name), value=VE.SP)
    assign(x=paste0("PVal_",name), value=PVal.SP)
    assign(x=paste0("Power_",name), value=Power.SP)
    assign(x=paste0("Res_",name), value=Res.SP)
  }
}

VE.Medians <- as.matrix(VE.Medians)
PVal.All <- as.matrix(PVal.All)
Power.All <- as.matrix(Power.All)

# for (scen in scens) {
#   for (param in params) {
#     name <- paste0("S",scen,"_P",param)
#     paramRow <- Varying[((scen-1)*12+as.numeric(param)),]
#     Power.SP <- get(paste0("Power_",name))
#     Res.SP <- get(paste0("Res_",name))
#     xlimval=c(0,100)
#     setEPS()
#     postscript(file=paste0(figfolder,"/Power/",name,".eps"),
#                width=6,height=6, paper="special")
#     plot(x=Power.SP[c("iRCT","cRCT_gamma_coxph","MEM","CPI",
#                              "MEM.Ex1","CPI.Ex1","lRR.NP.MeanPd","lRR.NP.ToCaWt",
#                              "lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1")]*100,
#          y=c(.9,.8,.7,.65,.55,.5,.4,.35,.25,.2), yaxt='n',ylab="",
#          type="p", pch=15,
#          ylim=c(.15,1),xlim=xlimval,xlab="Power (%)",
#          main=bquote("Power,"~beta==.(paramRow$Beta)*","~"VE"==.(paramRow$VE)))
#     for (j in seq(0,100,by=10)) {
#       abline(v=j, col="lightgray", lty="dotted")
#     }
#     if (paramRow$VE==0) {
#       abline(v=5, col="red", lty=4)
#     }
#     text(x=mean(xlimval), y=.95, labels="IRT")
#     text(x=mean(xlimval), y=.85, labels="CRT: Gamma Frailty")
#     text(x=mean(xlimval), y=.75, label="SW-CRT, MEM/CPI: Include First Period")
#     text(x=mean(xlimval), y=.6, label="SW-CRT, MEM/CPI: Exclude First Period")
#     text(x=mean(xlimval), y=.45, labels="SW-CRT, NPWP: Include First Period")
#     text(x=mean(xlimval), y=.3, labels="SW-CRT, NPWP: Exclude First Period")
#     dev.off()
#   }
# }

for (param in params) {
  name3 <- paste0("S",3,"_P",param)
  name4 <- paste0("S",4,"_P",param)
  paramRow3 <- Varying[((3-1)*12+as.numeric(param)),]
  paramRow4 <- Varying[((4-1)*12+as.numeric(param)),]
  Res3 <- get(paste0("Res_",name3))
  Res4 <- get(paste0("Res_",name4))
  
  if (paramRow3$VE==0) {
    xlimval=c(-.6,.6)
  } else {
    xlimval=c(0,1)
  }
  typenames <- c("iRCT","cRCT_gamma_coxph","MEM.Ex1","CPI.Ex1",
                 "lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
                 "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
                 "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")
  typelabs <- c("IRT","CRT","MEM","CPI","NPWP-1","NPWP-2","SC-1","SC-2","SC-Wt-1","SC-Wt-2")
  ys <- seq(from=1,by=-.05,length.out=length(typenames))
  setEPS()
  postscript(file=paste0(figfolder,"/Tx Estimates/Ests_P",param,".eps"),
             width=12,height=6, paper="special")
  par(mfrow=c(1,2))
  plot(x=Res3["50%",typenames],
       y=ys, xlim=xlimval, ylim=c(min(ys)-.05,1),
       type="p", pch=18, cex=1.8, yaxt='n', ylab="",
       xlab="Estimated VE",
       main=bquote("A)"~"VE Estimates,"~R[0]==.(round(paramRow3$R0,2))*","~"VE"==.(paramRow3$VE)*", SWT-A"))
  for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
    abline(v=j, col="lightgray", lty="dotted")
  }
  for (i in 1:length(typenames)) {
    abline(h=ys[i], col="lightgray", lty="dotted")
  }
  abline(v=paramRow3$VE, lty=5, col="gray", lwd=2)
  points(x=Res3["25%",typenames], y=ys, pch=20, cex=1.2)
  points(x=Res3["50%",typenames], y=ys, pch=18, cex=1.8)
  points(x=Res3["75%",typenames], y=ys, pch=20, cex=1.2)
  for (i in 1:length(typenames)) {
    lines(x=Res3[c("25%","75%"),typenames[i]], y=rep(ys[i],2), lty=1)
  }
  axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
  legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
         pt.cex=c(1.8,1.2), bg="white")
  
  plot(x=Res4["50%",typenames],
       y=ys, xlim=xlimval, ylim=c(min(ys)-.05,1),
       type="p", pch=18, cex=1.8, yaxt='n', ylab="",
       xlab="Estimated VE",
       main=bquote("B)"~"VE Estimates,"~R[0]==.(round(paramRow4$R0,2))*","~"VE"==.(paramRow4$VE)*", SWT-B"))
  for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
    abline(v=j, col="lightgray", lty="dotted")
  }
  for (i in 1:length(typenames)) {
    abline(h=ys[i], col="lightgray", lty="dotted")
  }
  abline(v=paramRow4$VE, lty=5, col="gray", lwd=2)
  points(x=Res4["25%",typenames], y=ys, pch=20, cex=1.2)
  points(x=Res4["50%",typenames], y=ys, pch=18, cex=1.8)
  points(x=Res4["75%",typenames], y=ys, pch=20, cex=1.2)
  for (i in 1:length(typenames)) {
    lines(x=Res4[c("25%","75%"),typenames[i]], y=rep(ys[i],2), lty=1)
  }
  axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
  legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
         pt.cex=c(1.8,1.2), bg="white")
  
  # 
  # plot(x=Res3[c("25%","50%","75%"),"iRCT"], y=rep(.95,3), xlim=xlimval, ylim=c(-.55,1),
  #      type="l", lty=1,
  #      xlab="Estimated VE", yaxt='n',ylab="",
  #      main=bquote("a."~"VE Estimates,"~R[0]==.(round(paramRow3$R0,2))*","~"VE"==.(paramRow3$VE)*", SWT-A"))
  # for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
  #   abline(v=j, col="lightgray", lty="dotted")
  # }
  # abline(v=paramRow3$VE, lty=2, col="red")
  # text(x=mean(xlimval), y=1, labels="IRT")
  # points(x=Res3[c("25%","50%","75%"),"iRCT"], y=rep(.95,3), pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.85, labels="CRT")
  # lines(x=Res3[c("25%","50%","75%"),"cRCT_gamma_coxph"], y=rep(.8,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"cRCT_gamma_coxph"], y=rep(.8,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.7, label="MEM")
  # lines(x=Res3[c("25%","50%","75%"),"MEM.Ex1"], y=rep(.65,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"MEM.Ex1"], y=rep(.65,3),
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval),y=.55, label="CPI")
  # lines(x=Res3[c("25%","50%","75%"),"CPI.Ex1"], y=rep(.5,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"CPI.Ex1"], y=rep(.5,3),
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.4, labels="NPWP-1")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.NP.MeanPd.Ex1"], y=rep(.35,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.NP.MeanPd.Ex1"], y=rep(.35,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.25, labels="NPWP-2")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.NP.ToCaWt.Ex1"], y=rep(.2,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.NP.ToCaWt.Ex1"], y=rep(.2,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.1, labels="SC-1")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.SC.MeanPd.Ex1"], y=rep(.05,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.SC.MeanPd.Ex1"], y=rep(.05,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=-.05, labels="SC-2")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.SC.ToCaWt.Ex1"], y=rep(-.1,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.SC.ToCaWt.Ex1"], y=rep(-.1,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=-.2, labels="SC-Wt-1")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.SCWt.MeanPd.Ex1"], y=rep(-.25,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.SCWt.MeanPd.Ex1"], y=rep(-.25,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=-.35, labels="SC-Wt-2")
  # lines(x=Res3[c("25%","50%","75%"),"lRR.SCWt.ToCaWt.Ex1"], y=rep(-.4,3), lty=1)
  # points(x=Res3[c("25%","50%","75%"),"lRR.SCWt.ToCaWt.Ex1"], y=rep(-.4,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
  #        pt.cex=c(1.5,1), bg="white")
  # 
  # plot(x=Res4[c("25%","50%","75%"),"iRCT"], y=rep(.9,3), xlim=xlimval, ylim=c(-.15,1),
  #      type="l", lty=1,
  #      xlab="Estimated VE", yaxt='n',ylab="",
  #      main=bquote("b."~"VE Estimates,"~R[0]==.(round(paramRow4$R0,2))*","~"VE"==.(paramRow4$VE)*", SWT-B"))
  # for (j in seq(from=xlimval[1],to=xlimval[2],by=.1)) {
  #   abline(v=j, col="lightgray", lty="dotted")
  # }
  # abline(v=paramRow4$VE, lty=2, col="red")
  # text(x=mean(xlimval), y=.95, labels="IRT")
  # points(x=Res4[c("25%","50%","75%"),"iRCT"], y=rep(.9,3), pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.85, labels="CRT")
  # lines(x=Res4[c("25%","50%","75%"),"cRCT_gamma_coxph"], y=rep(.8,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"cRCT_gamma_coxph"], y=rep(.8,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.75, label="MEM")
  # lines(x=Res4[c("25%","50%","75%"),"MEM.Ex1"], y=rep(.7,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"MEM.Ex1"], y=rep(.7,3),
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval),y=.65, label="CPI")
  # lines(x=Res4[c("25%","50%","75%"),"CPI.Ex1"], y=rep(.6,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"CPI.Ex1"], y=rep(.6,3),
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.55, labels="NPWP-1")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.NP.MeanPd.Ex1"], y=rep(.5,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.NP.MeanPd.Ex1"], y=rep(.5,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.45, labels="NPWP-2")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.NP.ToCaWt.Ex1"], y=rep(.4,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.NP.ToCaWt.Ex1"], y=rep(.4,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.35, labels="SC-1")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.SC.MeanPd.Ex1"], y=rep(.3,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.SC.MeanPd.Ex1"], y=rep(.3,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.25, labels="SC-2")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.SC.ToCaWt.Ex1"], y=rep(.2,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.SC.ToCaWt.Ex1"], y=rep(.2,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.15, labels="SC-Wt-1")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.SCWt.MeanPd.Ex1"], y=rep(.1,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.SCWt.MeanPd.Ex1"], y=rep(.1,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # text(x=mean(xlimval), y=.05, labels="SC-Wt-2")
  # lines(x=Res4[c("25%","50%","75%"),"lRR.SCWt.ToCaWt.Ex1"], y=rep(0,3), lty=1)
  # points(x=Res4[c("25%","50%","75%"),"lRR.SCWt.ToCaWt.Ex1"], y=rep(0,3), 
  #        pch=c(20,18,20), cex=c(1,1.5,1))
  # legend(x="bottom", legend=c("Median","Quartiles"), pch=c(18,20), horiz=TRUE, 
  #        pt.cex=c(1.5,1), bg="white")
  dev.off()
}

hm.colors <- snow.colors(n=72)[65:2]
hm.colors.gray <- gray.colors(n=64, start=.8, end=.2)
for (param in params) {
  name3 <- paste0("S",3,"_P",param)
  name4 <- paste0("S",4,"_P",param)
  paramRow3 <- Varying[((3-1)*12+as.numeric(param)),]
  paramRow4 <- Varying[((4-1)*12+as.numeric(param)),]
  fullResComp3 <- get(paste0("fullResComp_",name3))
  fullResComp4 <- get(paste0("fullResComp_",name4))
  
  xlimval3 <- c(2.5,10.5)
  xseq3 <- seq(.5,11.5,by=1)
  xlabels3 <- 3:10
  xlimval4 <- c(5.5,42.5)
  xseq4 <- seq(0.5,44.5,by=1)
  xlabels4 <- seq(6,42,by=2)
  Heat3 <- 1-exp(apply(fullResComp3[c("SCWt.lRR.Ex1","SC.lRR.Ex1","NP.lRR.Ex1"),,],
                      c(1,2), mean, na.rm=TRUE))
  Heat3 <- ifelse(Heat3 < -.2, -.2, Heat3)
  Heat4 <- 1-exp(apply(fullResComp4[c("SCWt.lRR.Ex1","SC.lRR.Ex1","NP.lRR.Ex1"),,],
                       c(1,2), mean, na.rm=TRUE))
  Heat4 <- ifelse(Heat4 < -.2, -.2, Heat4)
  setEPS()
  postscript(file=paste0(figfolder,"/Pd-Specific Heat Maps/Heat_P",param,".eps"),
             width=12,height=12, paper="special")
  par(mfrow=c(2,1))
  image.plot(x=xseq3,y=1:4,z=t(Heat3),xlim=xlimval3,ylim=c(1,4),
             ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
             main=bquote("A)"~"Estimated VE by Period, SWT-A"),
             col=hm.colors)
  axis(side=2,
       labels=c("SC-Wt","SC","NPWP"),
       at=c(1.5,2.5,3.5))
  axis(side=1, labels=xlabels3, at=xlabels3)
  abline(h=2, col=1, lty=1, lwd=1.5)
  abline(h=3, col=1, lty=1, lwd=1.5)
  
  image.plot(x=xseq4,y=1:4,z=t(Heat4),xlim=xlimval4,ylim=c(1,4),
             ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
             main=bquote("B)"~"Estimated VE by Period, SWT-B"),
             col=hm.colors)
  axis(side=2,
       labels=c("SC-Wt","SC","NPWP"),
       at=c(1.5,2.5,3.5))
  axis(side=1, labels=xlabels4, at=xlabels4)
  abline(h=2, col=1, lty=1, lwd=1.5)
  abline(h=3, col=1, lty=1, lwd=1.5)
  dev.off()
  
  setEPS()
  postscript(file=paste0(figfolder,"/Pd-Specific Heat Maps/GrayHeat_P",param,".eps"),
             width=12,height=12, paper="special")
  par(mfrow=c(2,1))
  image.plot(x=xseq3,y=1:4,z=t(Heat3),xlim=xlimval3,ylim=c(1,4),
             ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
             main=bquote("A)"~"Estimated VE by Period, SWT-A"),
             col=hm.colors.gray)
  axis(side=2,
       labels=c("SC-Wt","SC","NPWP"),
       at=c(1.5,2.5,3.5))
  axis(side=1, labels=xlabels3, at=xlabels3)
  abline(h=2, col=1, lty=1, lwd=1.5)
  abline(h=3, col=1, lty=1, lwd=1.5)
  
  image.plot(x=xseq4,y=1:4,z=t(Heat4),xlim=xlimval4,ylim=c(1,4),
             ylab="Analysis Method",yaxt='n',xaxt='n', xlab="Period",
             main=bquote("B)"~"Estimated VE by Period, SWT-B"),
             col=hm.colors.gray)
  axis(side=2,
       labels=c("SC-Wt","SC","NPWP"),
       at=c(1.5,2.5,3.5))
  axis(side=1, labels=xlabels4, at=xlabels4)
  abline(h=2, col=1, lty=1, lwd=1.5)
  abline(h=3, col=1, lty=1, lwd=1.5)
  dev.off()
}

for (param in c("02","03","05","06","08","09","11","12")) {
  paramRow3 <- Varying[((3-1)*12+as.numeric(param)),]
  beta <- paramRow3$Beta
  R0 <- paramRow3$R0
  paramRow4 <- Varying[((4-1)*12+as.numeric(param)),]
  
  S34Power <- Power.All[Power.All[,"Scen"] %in% c(3,4) & Power.All[,"Param"]==as.numeric(param),]
  S34TIE <- Power.All[Power.All[,"Scen"] %in% c(3,4) & Power.All[,"Beta"]==beta & Power.All[,"VE"]==0,]
  
  xlimval=c(0,100)
  typenames <- c("iRCT","cRCT_gamma_coxph","MEM.Ex1","CPI.Ex1",
                 "lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
                 "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
                 "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")
  typelabs <- c("IRT","CRT","MEM","CPI","NPWP-1","NPWP-2","SC-1","SC-2","SC-Wt-1","SC-Wt-2")
  ys <- seq(from=1,by=-.05,length.out=length(typenames))
  setEPS()
  postscript(file=paste0(figfolder,"/Power/Power_P",param,".eps"),
             width=12,height=6,paper="special")
  par(mfrow=c(1,2))
  plot(x=S34Power[S34Power[,"Scen"]==3,typenames]*100,y=ys,
       xlim=xlimval, ylim=c(min(ys)-.05,1),
       type="p", pch=c(15,16,rep(24,8)), bg=rep(1,10), 
       yaxt='n', ylab="",
       xlab="Power (%)",
       main=bquote("A)"~"Power,"~R[0]==.(round(R0,2))*","~"VE"==.(paramRow3$VE)))
  for (j in seq(from=xlimval[1],to=xlimval[2],by=10)) {
    abline(v=j, col="lightgray", lty="dotted")
  }
  for (i in 1:length(typenames)) {
    abline(h=ys[i], col="lightgray", lty="dotted")
  }
  points(x=S34Power[S34Power[,"Scen"]==3,typenames]*100, y=ys,
         pch=c(15,16,rep(24,8)), bg=rep(1,10), cex=1.5)
  points(x=S34Power[S34Power[,"Scen"]==4,typenames[3:6]]*100, y=ys[3:6],
         pch=25, bg=1, cex=1.5)
  if (param %in% c("02","03","05","06","08","09","11","12")) {
    points(x=S34Power[S34Power[,"Scen"]==4,typenames[7:10]]*100, y=ys[7:10],
           pch=25, bg=1, cex=1.5)
  }
  axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
  legend(x="bottom", legend=c("IRT","CRT","SWT-A","SWT-B"), pch=c(15,16,24,25), horiz=TRUE,
         pt.bg=rep(1,4), bg="white")
  
  plot(x=S34TIE[S34TIE[,"Scen"]==3,typenames]*100,y=ys,
       xlim=xlimval, ylim=c(min(ys)-.05,1),
       type="p", pch=c(15,16,rep(24,8)), bg=rep(1,10), 
       yaxt='n', ylab="",
       xlab="Type I Error (%)",
       main=bquote("B)"~"Type I Error,"~R[0]==.(round(R0,2))))
  for (j in seq(from=xlimval[1],to=xlimval[2],by=10)) {
    abline(v=j, col="lightgray", lty="dotted")
  }
  for (i in 1:length(typenames)) {
    abline(h=ys[i], col="lightgray", lty="dotted")
  }
  abline(v=5, lty=5, col="gray", lwd=2)
  points(x=S34TIE[S34TIE[,"Scen"]==3,typenames]*100, y=ys,
         pch=c(15,16,rep(24,8)), bg=rep(1,10), cex=1.5)
  points(x=S34TIE[S34TIE[,"Scen"]==4,typenames[3:6]]*100, y=ys[3:6],
         pch=25, bg=1, cex=1.5)
  if (param %in% c("02","03","05","06","08","09","11","12")) {
    points(x=S34TIE[S34TIE[,"Scen"]==4,typenames[7:10]]*100, y=ys[7:10],
           pch=25, bg=1, cex=1.5)
  }
  axis(side=2, at=ys, labels=typelabs, las=2, cex.axis=.8)
  legend(x="bottom", legend=c("IRT","CRT","SWT-A","SWT-B"), pch=c(15,16,24,25), horiz=TRUE,
         pt.bg=rep(1,4), bg="white")
  
  
  # plot(x=S34Power[S34Power[,"Scen"]==3,c("iRCT","cRCT_gamma_coxph",
  #                                        "MEM.Ex1","CPI.Ex1","lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
  #                                        "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
  #                                        "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")]*100,
  #      y=c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0), 
  #      main=bquote("a."~"Power,"~R[0]==.(round(R0,2))*","~"VE"==.(paramRow3$VE)),
  #      ylab="", yaxt="n", type="p",
  #      pch=c(15,16,rep(24,8)), bg=rep(1,10),
  #      ylim=c(-.15,1), xlim=xlimval, xlab="Power (%)")
  # points(x=S34Power[S34Power[,"Scen"]==4,c("MEM.Ex1","CPI.Ex1","lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1")]*100,
  #        y=c(.7,.6,.5,.4), pch=25, bg=1)
  # if (param %in% c("02","05","08","11")) {
  #   points(x=S34Power[S34Power[,"Scen"]==4,c("lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
  #                                            "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")]*100,
  #          y=c(.3,.2,.1,0), pch=25, bg=1)
  # }
  # for (j in c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0)) {
  #   abline(h=j, col=1, lty=1, lwd=.5)
  # }
  # for (j in seq(0,100,by=10)) {
  #   abline(v=j, col="lightgray", lty="dotted")
  # }
  # text(x=mean(xlimval), y=.95, labels="IRT")
  # text(x=mean(xlimval), y=.85, labels="CRT")
  # text(x=mean(xlimval), y=.75, label="MEM")
  # text(x=mean(xlimval), y=.65, label="CPI")
  # text(x=mean(xlimval), y=.55, labels="NPWP-1")
  # text(x=mean(xlimval), y=.45, labels="NPWP-2")
  # text(x=mean(xlimval), y=c(.35,.25,.15,.05), labels=c("SC-1","SC-2","SC-Wt-1","SC-Wt-2"))
  # legend(x="bottom", legend=c("IRT","CRT","SWT-A","SWT-B"), pch=c(15,16,24,25), horiz=TRUE,
  #        pt.bg=rep(1,4), bg="white")
  # 
  # plot(x=S34TIE[S34TIE[,"Scen"]==3,c("iRCT","cRCT_gamma_coxph",
  #                                    "MEM.Ex1","CPI.Ex1","lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1",
  #                                    "lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
  #                                    "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")]*100,
  #      y=c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0), yaxt='n',ylab="",
  #      type="p", pch=c(15,16,rep(24,8)), bg=rep(1,10),
  #      ylim=c(-.15,1),xlim=xlimval,xlab="Type I Error (%)",
  #      main=bquote("b."~"Type I Error,"~R[0]==.(round(R0,2))))
  # points(x=S34TIE[S34TIE[,"Scen"]==4,c("MEM.Ex1","CPI.Ex1","lRR.NP.MeanPd.Ex1","lRR.NP.ToCaWt.Ex1")]*100,
  #        y=c(.7,.6,.5,.4), pch=25, bg=1)
  # if (param %in% c("02","05","08","11")) {
  #   points(x=S34TIE[S34TIE[,"Scen"]==4,c("lRR.SC.MeanPd.Ex1","lRR.SC.ToCaWt.Ex1",
  #                                            "lRR.SCWt.MeanPd.Ex1","lRR.SCWt.ToCaWt.Ex1")]*100,
  #          y=c(.3,.2,.1,0), pch=25, bg=1)
  # }
  # for (j in c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0)) {
  #   abline(h=j, col=1, lty=1, lwd=.5)
  # }
  # for (j in seq(0,100,by=10)) {
  #   abline(v=j, col="lightgray", lty="dotted")
  # }
  # abline(v=5, col="red", lty=4)
  # text(x=mean(xlimval), y=.95, labels="IRT")
  # text(x=mean(xlimval), y=.85, labels="CRT")
  # text(x=mean(xlimval), y=.75, label="MEM")
  # text(x=mean(xlimval), y=.65, label="CPI")
  # text(x=mean(xlimval), y=.55, labels="NPWP-1")
  # text(x=mean(xlimval), y=.45, labels="NPWP-2")
  # text(x=mean(xlimval), y=c(.35,.25,.15,.05), labels=c("SC-1","SC-2","SC-Wt-1","SC-Wt-2"))
  # legend(x="bottom", legend=c("IRT","CRT","SWT-A","SWT-B"), pch=c(15,16,24,25), pt.bg=rep(1,4),
  #        horiz=TRUE, bg="white")
  dev.off()
}


for (VE in c(0,0.6,0.8)) {
  Medians.3 <- VE.Medians[VE.Medians[,"Scen"]==3 & VE.Medians[,"VE"]==VE,]
  Power.3 <- Power.All[Power.All[,"Scen"]==3 & Power.All[,"VE"]==VE,]
  Medians.4 <- VE.Medians[VE.Medians[,"Scen"]==4 & VE.Medians[,"VE"]==VE,]
  Power.4 <- Power.All[Power.All[,"Scen"]==4 & Power.All[,"VE"]==VE,]
  
  if (VE==0) {
    medianylim=c(-.5,.5)
    powerylim=c(0,100)
    medgrids=seq(-.4,.4,by=.2)
    powergrids=seq(0,100,by=20)
  } else if (VE==0.8) {
    medianylim=c(0,1)
    powerylim=c(0,100)
    medgrids=seq(0,1,by=.2)
    powergrids=seq(0,100,by=20)
  } else {
    medianylim=c(0,1)
    powerylim=c(0,100)
    medgrids=seq(0,1,by=.2)
    powergrids=seq(0,100,by=20)
  }
  setEPS()
  postscript(file=paste0(figfolder,"/Trends/Trend_VE",VE*100,".eps"),
      width=12,height=6,paper="special")
  par(mfrow=c(1,2))
  plot(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],ylim=medianylim,
       type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Median"~"VE Estimate"),
       main=bquote("A)"~"Median VE Estimate, VE"==.(VE)))
  points(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],pch=15)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,lty=2)
  points(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,pch=16)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.NP.ToCaWt.Ex1"],col=1,lty=3)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.NP.ToCaWt.Ex1"],col=1,pch=24,bg=1)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.NP.ToCaWt.Ex1"],col=1,lty=4)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.NP.ToCaWt.Ex1"],col=1,pch=25,bg=1)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SC.ToCaWt.Ex1"],col=1,lty=5)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SC.ToCaWt.Ex1"],col=1,pch=2)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SC.ToCaWt.Ex1"],col=1,lty=6)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SC.ToCaWt.Ex1"],col=1,pch=6)
  for (j in medgrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  abline(h=VE, lty=5, col="gray", lwd=2)
  legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A NPWP-2", "SWT-A SC-2","CRT","SWT-B NPWP-2","SWT-B SC-2"),
         lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  
  plot(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,ylim=powerylim,type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Power"~"(%)"),
       main=bquote("B)"~"Power, VE"==.(VE)))
  points(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,pch=15)
  lines(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,lty=2)
  points(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,pch=16, cex=1.3)
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,lty=3)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,pch=24, bg=1)
  lines(x=Power.4[,"R0"],y=Power.4[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,lty=4)
  points(x=Power.4[,"R0"],y=Power.4[,"lRR.NP.ToCaWt.Ex1"]*100,col=1,pch=25, bg=1)
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,lty=5)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,pch=2)
  if (VE==0.6 | VE==0.8) {
    lines(x=Power.4[,"R0"],y=Power.4[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,lty=6)
    points(x=Power.4[,"R0"],y=Power.4[,"lRR.SC.ToCaWt.Ex1"]*100,col=1,pch=6)
  }
  for (j in powergrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  if (VE==0) {
    abline(h=5, lty=5, col="gray", lwd=2)
    legend(x="bottomleft",ncol=2,legend=c("IRT","SWT-A NPWP-2", "SWT-A SC-2","CRT","SWT-B NPWP-2"),
           lty=c(1,3,5,2,4),pch=c(15,24,2,16,25),pt.bg=c(NA,1,NA,NA,1), bg="white")
  } else {
    legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A NPWP-2", "SWT-A SC-2","CRT","SWT-B NPWP-2","SWT-B SC-2"),
           lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  }
  dev.off()
  
  setEPS()
  postscript(file=paste0(figfolder,"/Trends/AltTrend1_VE",VE*100,".eps"),
             width=12,height=6,paper="special")
  par(mfrow=c(1,2))
  plot(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],ylim=medianylim,
       type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Median"~"VE Estimate"),
       main=bquote("A)"~"Median VE Estimate, VE"==.(VE)))
  points(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],pch=15)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,lty=2)
  points(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,pch=16)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.NP.MeanPd.Ex1"],col=1,lty=3)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.NP.MeanPd.Ex1"],col=1,pch=24,bg=1)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.NP.MeanPd.Ex1"],col=1,lty=4)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.NP.MeanPd.Ex1"],col=1,pch=25,bg=1)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SC.MeanPd.Ex1"],col=1,lty=5)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SC.MeanPd.Ex1"],col=1,pch=2)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SC.MeanPd.Ex1"],col=1,lty=6)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SC.MeanPd.Ex1"],col=1,pch=6)
  for (j in medgrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  abline(h=VE, lty=5, col="gray", lwd=2)
  legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A NPWP-1", "SWT-A SC-1","CRT","SWT-B NPWP-1","SWT-B SC-1"),
         lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  
  plot(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,ylim=powerylim,type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Power"~"(%)"),
       main=bquote("B)"~"Power, VE"==.(VE)))
  points(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,pch=15)
  lines(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,lty=2)
  points(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,pch=16, cex=1.3)
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.NP.MeanPd.Ex1"]*100,col=1,lty=3)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.NP.MeanPd.Ex1"]*100,col=1,pch=24, bg=1)
  lines(x=Power.4[,"R0"],y=Power.4[,"lRR.NP.MeanPd.Ex1"]*100,col=1,lty=4)
  points(x=Power.4[,"R0"],y=Power.4[,"lRR.NP.MeanPd.Ex1"]*100,col=1,pch=25, bg=1)
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SC.MeanPd.Ex1"]*100,col=1,lty=5)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.SC.MeanPd.Ex1"],col=1,pch=2)
  if (VE==0.6 | VE==0.8) {
    lines(x=Power.4[,"R0"],y=Power.4[,"lRR.SC.MeanPd.Ex1"]*100,col=1,lty=6)
    points(x=Power.4[,"R0"],y=Power.4[,"lRR.SC.MeanPd.Ex1"],col=1,pch=6)
  }
  for (j in powergrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  if (VE==0) {
    abline(h=5, lty=5, col="gray", lwd=2)
    legend(x="bottomleft",ncol=2,legend=c("IRT","SWT-A NPWP-1", "SWT-A SC-1","CRT","SWT-B NPWP-1","SWT-B SC-1"),
           lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  } else {
    legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A NPWP-1", "SWT-A SC-1","CRT","SWT-B NPWP-1","SWT-B SC-1"),
           lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  }
  dev.off()
  
  setEPS()
  postscript(file=paste0(figfolder,"/Trends/AltTrendSCWt_VE",VE*100,".eps"),
             width=12,height=6,paper="special")
  par(mfrow=c(1,2))
  plot(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],ylim=medianylim,
       type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Median"~"VE Estimate"),
       main=bquote("A)"~"Median VE Estimate, VE"==.(VE)))
  points(x=Medians.3[,"R0"],y=Medians.3[,"iRCT"],pch=15)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,lty=2)
  points(x=Medians.3[,"R0"],y=Medians.3[,"cRCT_gamma_coxph"],col=1,pch=16)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SCWt.MeanPd.Ex1"],col=1,lty=3)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SCWt.MeanPd.Ex1"],col=1,pch=24,bg=1)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SCWt.MeanPd.Ex1"],col=1,lty=4)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SCWt.MeanPd.Ex1"],col=1,pch=25,bg=1)
  lines(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,lty=5)
  points(x=Medians.3[,"R0"],y=Medians.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,pch=2)
  lines(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SCWt.ToCaWt.Ex1"],col=1,lty=6)
  points(x=Medians.4[,"R0"],y=Medians.4[,"lRR.SCWt.ToCaWt.Ex1"],col=1,pch=6)
  for (j in medgrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  abline(h=VE, lty=5, col="gray", lwd=2)
  legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A SC-Wt-1", "SWT-A SC-Wt-2","CRT","SWT-B SC-Wt-1","SWT-B SC-Wt-2"),
         lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  
  plot(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,ylim=powerylim,type="l",col=1,lty=1,
       xlab=bquote(R[0]),ylab=bquote("Power"~"(%)"),
       main=bquote("B)"~"Power, VE"==.(VE)))
  points(x=Power.3[,"R0"],y=Power.3[,"iRCT"]*100,pch=15)
  lines(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,lty=2)
  points(x=Power.3[,"R0"],y=Power.3[,"cRCT_gamma_coxph"]*100,col=1,pch=16, cex=1.3)
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.MeanPd.Ex1"]*100,col=1,lty=3)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.MeanPd.Ex1"]*100,col=1,pch=24, bg=1)
  if (VE==0.6 | VE==0.8) {
    lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.MeanPd.Ex1"]*100,col=1,lty=4)
    points(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.MeanPd.Ex1"]*100,col=1,pch=25, bg=1)
  }
  lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,lty=5)
  points(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,pch=2)
  if (VE==0.6 | VE==0.8) {
    lines(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,lty=6)
    points(x=Power.3[,"R0"],y=Power.3[,"lRR.SCWt.ToCaWt.Ex1"],col=1,pch=6)
  }
  for (j in powergrids) {
    abline(h=j, col="lightgray", lty="dotted")
  }
  if (VE==0) {
    abline(h=5, lty=5, col="gray", lwd=2)
    legend(x="bottomleft",ncol=2,legend=c("IRT","SWT-A SC-Wt-1", "SWT-A SC-Wt-2","CRT","SWT-B SC-Wt-1","SWT-B SC-Wt-2"),
           lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  }
  else {
    legend(x="bottomright",ncol=2,legend=c("IRT","SWT-A SC-Wt-1", "SWT-A SC-Wt-2","CRT","SWT-B SC-Wt-1","SWT-B SC-Wt-2"),
           lty=c(1,3,5,2,4,6),pch=c(15,24,2,16,25,6),pt.bg=c(NA,1,NA,NA,1,NA), bg="white")
  }
  dev.off()
}





