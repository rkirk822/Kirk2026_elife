outlist<-list()
folders <-
  list.files("Documents/EPhysTransfer/Klf913_FI/",
             pattern = "Cell[0-9]+$",
             full.names = T)
for (f in folders) {
  inFold <- f
  print(inFold)
  flist <- list.files(path=inFold, recursive = F, pattern=extPatt)
  if (length(flist) > 0) {    
    bMultFIs <- F
    folds <- inFold
    print(paste("single FI mode"))
  } else {
    bMultFIs <- T             # no data files, so maybe this is a folder of folders
    folds <- list.dirs(path=inFold, recursive = F)
    print(paste("multi FI mode"))
  }
  for (fold in folds) {
    print(paste("folder: ",fold))
    flist <- list.files(path=fold, recursive = F, pattern=extPatt)
    iFold <- which(folds==fold)
    fName <- flist[1]
    #
    #   Create folder for results
    #
    outdir1 <- paste0(fold,"/results_",Sys.Date(),"/")
    bExit=F
    resNum = 0
    while (bExit==F) {
      if (dir.exists(outdir1)==FALSE) {
        outdir <- outdir1
        dir.create(outdir)
        bExit=T
      } else {
        resNum = resNum+1
        outdir1 = paste0(fold,"/results_",Sys.Date(),letters[resNum],"/")
        if (resNum >25) {
          stop("too many results files")
        }
      }
    }
    #
    # test first file name for 4 digit format for cells, trials
    #
    if(substr(fName,1,5)=="Cell_" & (substr(fName,10,10)=="_") & (substr(fName,15,15)==".")==F) {
      msg <- paste("filename format error: ",fName)
      abort(msg)
    }
    nFiles <- length(flist) 
    if (nFiles < 1) stop("no data files")
    
    cell <- vector(mode = "integer", length = nFiles)
    trial <- vector(mode = "integer", length = nFiles)
    condNum <- vector(mode = "integer", length = nFiles)  # this is an integer index into steps, or anything else varied across trials
    titl <- vector(mode = "character", length = nFiles)
    rin <- vector(mode = "numeric", length = nFiles)
    tau <- vector(mode = "numeric", length = nFiles)
    vRest <- vector(mode = "numeric", length = nFiles)
    vFinal <- vector(mode = "numeric", length = nFiles)
    stepAmp <- vector(mode = "numeric", length = nFiles)
    nSpikes <- vector(mode = "integer", length = nFiles)
    rateMn <- vector(mode = "numeric", length = nFiles)
    rate12 <- vector(mode = "numeric", length = nFiles)
    rate23 <- vector(mode = "numeric", length = nFiles)
    fPath <- vector(mode = "character", length = nFiles)
    
    for (iFile in 1:nFiles) {
      fName <- flist[iFile]
      fPath[iFile] <- paste0(fold,"/",fName)
      cell[iFile] <- as.integer(substr(fName,6,9))
      trial[iFile] <- as.integer(substr(fName,11,14))
      titl[iFile] <- paste0("C",cell[iFile],"T",trial[iFile])
      print(paste(iFile,":",titl[iFile]))
      dList <- getIgorBin(fPath[iFile])
      y <- dList$data
      dkeyVals <- dList$keyVals
      if (mean(y[1:100]) > -0.01){ # cludge to correct files that are off by factor of 10 (e.g. rest is -7mv instead of -70)
        y = y*10    
      }
      
      #
      # get the step parameters and step amplitude from each file
      #
      if ("SealTestAmp(V)" %in% dkeyVals$keys) {
        print(paste("voltage clamp file. Zero spikes recorded for this trial."))
        next
      }
      stepParams <- getStepParams(dkeyVals)
      condNum[iFile] <- getStimNum(stepParams)
      stepAmp[iFile] <- as.numeric(stepParams[keys=="stAmp", vals])
      stepDur <- as.numeric(stepParams[keys=="stWidth", vals])
      passiveProps <- sealRinIclmp(y,dkeyVals, defaultSR)
      rin[iFile] <- passiveProps[keys=="rIn", vals]
      vRest[iFile] <- passiveProps[keys=="VmInit", vals]
      vFinal[iFile] <- passiveProps[keys=="VmFinal", vals]
      rin[iFile] <- passiveProps[keys=="rIn", vals]
      rin[iFile] <- passiveProps[keys=="rIn", vals]
      
      dur <- as.numeric(dkeyVals[keys=="Length(s)", vals])
      sr <- 1/as.numeric(dkeyVals[keys=="XDelta(s)", vals])
      trace <- data.table(t=c(1:(dur*sr)), y=y)
      
      riseLocs <- detectSpikes2(y, akeyVals)
      if (bDisplay) {
        fig <- plot_ly(trace, x=~t, y=~y, name='trace', type='scatter', mode='lines')
        fig <- fig %>% add_trace(x=~riseLocs, y=~y[riseLocs], name='rises', mode='markers')
        l <- layout(fig, title = titl[iFile])
        print(l)
      }
      
      
      if(sum(riseLocs)==0) {   # if there are no spikes, riseLocs set to scalar zero in detectSpikes
        if (bVerbose) print(paste("no spikes"))
      } else {
        spikeList <- measureSpikes(y, riseLocs, dkeyVals, akeyVals, stepParams)
        spikeProps <- spikeList[[1]]
        if (bVerbose) print(paste(nrow(spikeProps)," spikes"))
        outf <- paste0(outdir,titl[iFile],"_spikeProps.txt")
        fwrite(spikeProps,outf)
        
        #
        # spikeProps has one row per spike, cols are: pkLocs, pkAmps, minLocs, minAmps, pkRiseLocs, PkRiseAmps, threshLocs, threshAmps,ISI
        # ISI is sec, -Locs are in samples, Amps are in V.
        #
        nSpikes[iFile] <- nrow(spikeProps)
        nGoodSpikes <- sum(spikeProps$ISI > 0) # # of spikes with adequate repolarization
        rateMn[iFile] <- nGoodSpikes / stepDur
        if (nGoodSpikes >1) {
          rate12[iFile] <- round((1/spikeProps$ISI[2]),3)
        }
        if (nSpikes[iFile] >2) {
          rate23[iFile] <- round((1/spikeProps$ISI[3]),3)
        }
        spikeKeys <- spikeList[[2]]
        spikeVals <- spikeList[[3]]
      } # end if there are no spikes
    }
    
    resDT <- data.table(cell,condNum,titl,rin,tau,vRest,stepAmp,nSpikes,rateMn,rate12,rate23,fPath)
    #resDT<-mutate(resDT,Resting=ifelse(vRest>-0.072,"Normal (~-70mV)","Hyperpolarized (~-80mV)"))
    cellid<-basename(inFold)
    outlist[[cellid]]<-resDT

  }
}
#
# resDT is a data.table that has a calculated mean and inst. rate for each file; FIcurve is a summary of this that  contains means and SEM across repetitions
#
outdf<-do.call(rbind,outlist)


#For repairing miscoded cell in K9/13 data (P23) and getting rid of hyperpolarized trials
outdf$cell[which(grepl(outdf$titl,pattern = "C129T*")&as.numeric(str_split_fixed(outdf$titl,pattern = "T",2)[,2])>57)]<-130
outdf<-outdf%>%filter(grepl(Resting,pattern="Normal"))


#For eliminating botched trials in K9/13 data (P18)
outdf<-outdf%>%filter(!grepl(titl,pattern=paste0("C280T",as.character(seq(18)),collapse = "$|")))%>%
  filter(!grepl(titl,pattern=paste0("C282T",as.character(seq(16)),collapse = "$|")))


#Eliminate 'bad' cells, 
klfgood<-c(126:134,263,265,266,280,281,282,283,286,288)
klfpgood<-c(klfgood,259:261,264,284,285,287)
scrgood<-c(117:121,123,153,155,255,257,267,268,270:273)
scrpgood<-c(scrgood,269,274)

#outdf2<-outdf%>%filter(cell%in%c(klfpgood,scrpgood))%>%filter(-0.073<=vRest&vRest<=-0.068)
outdf2<-outdf%>%filter(Resting=="Normal (~-70mV)")%>%filter(-0.073<=vRest&vRest<=-0.068)
resDT<-outdf2%>%group_by(cell,stepAmp,Virus)%>%summarise(rateMn = mean(rateMn), rate12 = mean(rate12), rate23 = mean(rate23))

FIcurve <- resDT %>% group_by(stepAmp,Virus)%>%summarise(meanMn = mean(rateMn), mean12 = mean(rate12), mean23 = mean(rate23), 
                                                             n=n(), semMn=(sd(rateMn)/sqrt(n)), sem12=(sd(rate12)/sqrt(n)), sem23=(sd(rate23)/sqrt(n)), stepAmp = mean(stepAmp))

FIcurve%>%mutate(Virus=case_when(Virus=="Scr-mCherry|Scramble" ~ "Scramble gRNA", TRUE ~ Virus))%>%filter(n>3)%>%ggplot(aes(x=stepAmp*1e12,y=meanMn,color=Virus))+geom_line()+  
  geom_errorbar(aes(ymin=meanMn-semMn, ymax=meanMn+semMn), width=.00000000001)+
  theme_classic()+labs(y="Mean Firing Rate (spikes/s)",x="Current Amplitude (pA)",color="gRNA Target")+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, NA)) + scale_colour_manual(values = c("cyan","blue4"))+
  theme(text=element_text(size=15))

`#Combine P23 & P18 sets
outdf<-rbind(outdf23,outdf18)

outdf2<-outdf%>%filter(cell%in%c(klfpgood,scrpgood))%>%filter(-0.073<=vRest&vRest<=-0.068)
resDT<-outdf2%>%group_by(cell,stepAmp,Virus,Age)%>%summarise(rateMn = mean(rateMn), rate12 = mean(rate12), rate23 = mean(rate23))
FIcurve <- resDT %>% group_by(stepAmp,Virus,Age)%>%summarise(meanMn = mean(rateMn), mean12 = mean(rate12), mean23 = mean(rate23), 
                               n=n(), semMn=(sd(rateMn)/sqrt(n)), sem12=(sd(rate12)/sqrt(n)), sem23=(sd(rate23)/sqrt(n)), stepAmp = mean(stepAmp))

FIcurve%>%filter(n>3)%>%ggplot(aes(x=stepAmp*1e12,y=meanMn,color=Virus,linetype=Age))+geom_line()+  
  geom_errorbar(aes(ymin=meanMn-semMn, ymax=meanMn+semMn), width=.00000000001)+
  theme_classic()+labs(y="Mean Firing Rate (spikes/s)",x="Current Amplitude (pA)")+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, NA)) + scale_colour_manual(values = c("cadetblue3","blue4","red"))

outf <- paste0(outdir,titl[iFile],"_spikeRates.txt")
fwrite(resDT, outf)
outf <- paste0(outdir,titl[iFile],"_FIcurve.txt")
fwrite(FIcurve, outf)

p <- ggplot(resDT, aes(x=stepAmp, y=rateMn)) + geom_point() + geom_line(data=FIcurve, aes(x=stepAmp, y=meanMn))

outf <- paste0(outdir,titl[iFile],"_meanFI.pdf")

dev.new(pdf)
pdf(outf, width=10, height=6, pointsize=8)
plot(p)
dev.off()

q <- ggplot(resDT, aes(x=stepAmp, y=rate12))+geom_point()+geom_point(aes(y=rate23), color="red")+geom_line(data=FIcurve, aes(x=stepAmp, y=mean12))+geom_line(data=FIcurve, aes(x=stepAmp, y=mean23), color="red")

outf <- paste0(outdir,titl[iFile],"_instFI.pdf")

dev.new(pdf)
pdf(outf, width=10, height=6, pointsize=8)
plot(q)
dev.off()