
basefolder <- "/home/lk143/TimeVarID/OB_Sims"
res2folder <- paste0(basefolder,"/FullRes2")
NPSCfolder <- paste0(basefolder,"/NPSCRes")
Compfolder <- paste0(basefolder,"/Compiled")
load(paste0(basefolder,"/Varying.Rda"))


######################################################
###  function to easily use R packages on cluster  ###
######################################################
package.install = function(pack) {
  local({r <- getOption("repos");r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})
  
  # name of package to install / load
  pack = pack
  
  if (pack %in% rownames(installed.packages())) {
    library(pack, character.only=T)
  } else {
    if (pack %in% rownames(installed.packages(lib.loc='/home/lk143/apps/R_3.5.1/library'))) {
      library(pack, lib.loc='/home/lk143/apps/R_3.5.1/library', character.only=T)
    } else {
      install.packages(pack, lib='/home/lk143/apps/R_3.5.1/library')
      library(pack, lib.loc='/home/lk143/apps/R_3.5.1/library', character.only=T)
    }
  }
}


### Packages Needed: ###
library(abind)

scen <- as.numeric(substring(Sys.getenv('SLURM_JOB_NAME'), first=12, last=12))
for (scens in c(1,2,3,4)) {
  for (param in c("01","02","03","04","05","06","07","08","09","10","11","12")) {
    name <- paste0("S",scen,"_P",param)
    print(paste0("Beginning: ",name))
    TE.RD.Comp <- NULL
    TE.RR.Comp <- NULL
    TE.lRR.Comp <- NULL
    ResArray.Comp <- NULL
    errlist <- NULL
    for (num in 1:500) {
      name2 <- paste0(name,"_",num)
      if (file.exists(paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))) {
        load(paste0(res2folder,"/FullRes2_",name,"_",num,".Rda"))
        ResOut <- fullRes$CombRes
        ResArray.Comp <- abind(ResArray.Comp, ResOut, along=3)
        for (conttype in c("RD","RR","lRR")) {
          matname <- paste0("TE.",conttype,".num")
          assign(matname, c("Sim"=num))
          for (antype in c("NP","SC","SCWt")) {
            for (ex1 in c("","Ex1")) {
              row0ex1 <- ifelse(ex1=="Ex1",".Ex1","")
              row0 <- paste0(antype,".",conttype,row0ex1)
              innames <- names(get(matname))
              dat <- c(get(matname),
                       mean(ResOut[row0,], na.rm=TRUE),
                       sum(ResOut[row0,]*ResOut["IntClusts",], 
                           na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut[paste0("IntClusts",ex1),], 
                                           na.rm=TRUE),
                       sum(ResOut[row0,]*ResOut["ContCases",], 
                           na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut["ContCases",], 
                                           na.rm=TRUE),
                       sum(ResOut[row0,]*ResOut["TotCases",], 
                           na.rm=TRUE)/sum(ifelse(is.finite(ResOut[row0,]),1,0)*ResOut[paste0("TotCases",ex1),], 
                                           na.rm=TRUE))
              names(dat) <- c(innames, paste0(antype,".",c("MeanPd","MeanCP","CoCaWt","ToCaWt"),row0ex1))
              assign(matname, dat)
            }
          }
          assign(paste0("TE.",conttype,".Comp"), 
                 rbind(get(paste0("TE.",conttype,".Comp")), get(matname)))
        }
      } else {
        errlist <- c(errlist,num)
      }
    }
    print(paste0("Errors for sim:",name))
    print(paste(errlist,collapse=","))
    TE.RD <- TE.RD.Comp
    TE.RR <- TE.RR.Comp
    TE.lRR <- TE.lRR.Comp
    fullResComp <- ResArray.Comp
    save(TE.RD, file=paste0(Compfolder,"/TERD_",name,".Rda"))
    save(TE.RR, file=paste0(Compfolder,"/TERR_",name,".Rda"))
    save(TE.lRR, file=paste0(Compfolder,"/TElRR_",name,".Rda"))
    save(fullResComp, file=paste0(Compfolder,"/FullRes2Comp_",name,".Rda"))
  }
}