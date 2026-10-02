#
# ephysUtils.R  some utilities for analysis of IGOR electrophysiology data. 
# to use these in your code use: source("ephysUtils.R") if this file is in the same directory as your code, or
# source("directory_where_you_put_it/ephysUtils.R") if you put it in "/directory_where_you_put_it"
#

#************************************************************************
fileDisplay <- function(inFile="") {
  box::use(IgorR)
  if (inFile=="") {
    inFile <- rstudioapi::selectFile(
      caption = "Select Igor File",
      label = "Select",
      existing = TRUE)
  }
  print(paste(inFile))
  y <- read.ibw(inFile)

  t <- 1:(length(y))
  dat <- data.table(x=t, y=y)
  fig <- plot_ly(dat, x=~t, y=~y, name='orig', type='scatter', mode='lines') 
  
  print (fig)
  #browser()
  
}
#************************************************************************
fileAndDerivDisplay <- function(inFile="") {
  box::use(IgorR)
  if (inFile=="") {
    inFile <- rstudioapi::selectFile(
      caption = "Select Igor File",
      label = "Select",
      existing = TRUE)
  }
  print(paste(inFile))
  y <- read.ibw(inFile)
  dy <- y - shift(y, fill=y[1])               # 1st deriv
  
  t <- 1:(length(y))
  dat <- data.table(x=t, y=y, dy=dy)
  fig <- plot_ly(dat, x=~t, y=~y, name='orig', type='scatter', mode='lines') %>%
    add_trace(type="scatter", name="deriv", x=~t, y = ~dy*20, mode = 'lines') 
  
  print (fig)
  #browser()
  
}
#************************************************************************
derivDisplay <- function(inFile="",npts=1) {
  box::use(IgorR)
  box::use(RcppRoll)
  if (inFile=="") {
    inFile <- rstudioapi::selectFile(
      caption = "Select Igor File",
      label = "Select",
      existing = TRUE)
  }
  print(paste(inFile))
  y <- read.ibw(inFile)
  dy <- y - shift(y, fill=y[1])               # 1st deriv
  dyn <- roll_sum(dy,n=npts)
  
  t <- 1:(length(dyn))
  dat <- data.table(x=t,dyn=dyn)
  fig <- plot_ly(dat, x=~t, y=~dyn, type='scatter', mode='lines') 
  print (fig)
  #browser()
  
}#************************************************************************
histDisplay <- function(inFile="", start=1, end=-1) {
  box::use(IgorR)
  if (inFile=="") {
    inFile <- rstudioapi::selectFile(
      caption = "Select Igor File",
      label = "Select",
      existing = TRUE)
  }
  print(paste(inFile))
  y <- read.ibw(inFile)
  if (end==-1) end=length(y)
  ysub <- y[start:end]
  print(paste("start=",start,"end=",end))
  fig <- plot_ly(x =~ ysub, type="histogram")
#  dat <- data.table(x=t, y=y, dy=dy)
#  fig <- plot_ly(dat, x=~t, y=~y, name='orig', type='scatter', mode='lines') %>%
#    add_trace(type="scatter", name="deriv", x=~t, y = ~dy*20-0.07, mode = 'lines')
  
  print (fig)
  #browser()
  
}
#************************************************************************

sealRin <- function(xf, keyVals, defaultSR) {
  #
  # measures input resistance and holding current for v-clamp traces
  #
  sealStrt <- as.numeric(keyVals[keys=="SealTestStartX(s)", vals])
  sealLen <- as.numeric(keyVals[keys=="SealTestLength(s)", vals])
  sealAmp <- as.numeric(keyVals[keys=="SealTestAmp(V)", vals])
  #print(paste("sealAmp=",sealAmp))
  xDel <- as.numeric(keyVals[keys=="XDelta(s)", vals])            # = 1/sampling rate
  xLen <- as.numeric(keyVals[keys=="Length(s)", vals])      # inc seal test
  filesr <- 1/xDel
  if (filesr!=defaultSR) {
    stop(paste("SR for file is",filesr,"but default is ",defaultSR,": change sr default in main loop"))  
  }
  
  asymWin = 1 + sr*c((sealStrt+sealLen-0.1),(sealStrt+sealLen))
  sealWin = 1 + sr*c((sealStrt),(sealStrt+sealLen))
  baseEnd = sr*sealStrt
  bas = mean(xf[1:(baseEnd)])        # report in pA
  basPA = round(bas*1e12)        # report in pA
  asymp <- mean(xf[(asymWin[1]):(asymWin[2])])
  rin <- sealAmp/(asymp-bas)        # R = V/I; mV / pA = GOhms (but in Ohms)
  rinMO <-round(rin/1e6)            # report in megOhms
#  print(paste("asymp=",round(asymp*1e12),"from:",asymWin[1], "to", asymWin[2], " rin=",rinMO, "basPA=",basPA))
  return(c(basPA,rinMO))
}
#************************************************************************

sealRinIclmp <- function(x, keyVals, sr) {
  #
  # measures input resistance and resting Vm for I-clamp traces
  #
#browser()
  sealStrt <- as.numeric(keyVals[keys=="SealTestStartX(s)", vals])
  sealLen <- as.numeric(keyVals[keys=="SealTestLength(s)", vals])
  sealAmp <- as.numeric(keyVals[keys=="SealTestAmp(A)", vals])
  if((sealStrt==0) | (sealLen==0) | (sealAmp==0)) {
    warning(paste("sealAmp=",sealAmp,"strt=",sealStrt,"len=",sealLen))
    return()
  }
  xDel <- as.numeric(keyVals[keys=="XDelta(s)", vals])            # = 1/sampling rate
  xLen <- as.numeric(keyVals[keys=="Length(s)", vals])      # inc seal test
  filesr <- 1/xDel
  if (filesr!=sr) {
    stop(paste("SR for file is",filesr,"but default is ",sr,": change sr default in main loop"))  
  }
  
  asymWin = 1 + sr*c((sealStrt+sealLen-0.11),(sealStrt+sealLen-0.01))
#  sealWin = 1 + sr*c((sealStrt),(sealStrt+sealLen))
  baseEnd = (sr*sealStrt)-1
  bas = mean(x[1:(baseEnd)])        # should be in V
basPA = round(bas*1e12)        # 
  asymp <- mean(x[(asymWin[1]):(asymWin[2])])
  rin <- (asymp-bas)/sealAmp        # R = V/I; mV / pA = GOhms (but in Ohms)
  rinMO <-round((rin/1e6),3)            # report in megOhms
#print(paste("asymp=",round(asymp,3),"from:",asymWin[1], "to", asymWin[2], " rin=",rinMO, "bas=",bas))
  
###TODO fit exponential to get tau; get final Vm from end of trial
  keys <- c("VmInit", "VmFinal", "rIn", "tau")
  vals <- c(bas,0,rin,0)
  passiveProps <- data.table(keys,vals)
}
#************************************************************************
getIgorBin <- function(path) {
#
# Reads an Igor binary wave file and splits the wave note into a data.table consisting of a vector of keys 
# and vector of values, returns a list consisting of the keyVals table and the data
#
box::use(IgorR)
x <- read.ibw(path)
# get kw:value pairs out of wave note, which is an attr of the vector x called "Note".
notePrs <- str_split(attr(x, "Note"), ";", simplify=T)     # simplify returns char vec/matrix instead of list
pairNums <- c(1:length(notePrs))
splitPairs <- str_split(notePrs, ":", n=2, simplify = T)
paramNames <- splitPairs[pairNums]
paramVals <- splitPairs[(pairNums + length(notePrs))]
keyVals <- data.table(keys=paramNames, vals=paramVals)  # NOTE; these are strings, need to be converted

return (list(keyVals=keyVals, data=x))
}

#************************************************************************
#
# from wave note parsed into keys, vals
#
getStepParams <- function(keyVals) {
#browser()

  keys[1] <- "protocol"; vals[1] <- keyVals[keys=="StimProtocol", vals]
  keys[2] <- "stStart"; vals[2] <- keyVals[keys=="StepStart(s)", vals]
  keys[3] <- "stWidth"; vals[3]  <- keyVals[keys=="Width", vals]
  keys[4] <- "nStim"; vals[4]  <- keyVals[keys=="NStims", vals]
  keys[5] <- "st1Amp"; vals[5]  <- keyVals[keys=="1st Stim Amp.", vals]
  keys[6] <- "stAmpDiff"; vals[6]  <- keyVals[keys=="Amp. Diff", vals]
  keys[7] <- "stAmp"; vals[7]  <- keyVals[keys=="Stim Amp.", vals]
  
  stepParams <- data.table(keys,vals)
  
  return(stepParams)
}
  #************************************************************************
  #
  # Call after parsing wavenote into stepParams via getStepParams
  #
  getStimNum <- function(stepParams) {
  
  protocol <- stepParams[keys=="protocol", vals]
  stAmp <- as.numeric(stepParams[keys=="stAmp", vals])
  st1Amp <- as.numeric(stepParams[keys=="st1Amp", vals])
  stAmpDiff <- as.numeric(stepParams[keys=="stAmpDiff", vals])
  if (protocol=="NStimsAtFixedAmpDiff") {
    stimNum <- (stAmp - st1Amp)/stAmpDiff
  } else {
    if (protocol=="") {
      stimNum <- 1
    } else {
      print(paste("irregular protocol: ",protocol))
    } # endif single trial protocol
  } # endif protocol is NStimsAtFixedAmpDiff
  return(stimNum)                                           # 
}