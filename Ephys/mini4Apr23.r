###################################################################
#
# April 23rd: adding code to write detection/extraction parameters so that we can keep track of these and so they can be used
# by avgMini. Also recently split the summary results into two files (each for good and bad minis). The first file retains the same name:
# Cell#filtSumG.txt and Cell#filtSumB.txt. The second files: Cell#filtSumG2.txt and Cell#filtSumB2.txt include most of the p-values and
# residuals for the fits to baseline, rise and decay.
#
# April 5th
# changing extracted mini code so that minis are aligned on mid-rise at the time they are saved. They will also be baseline
# subtracted. This will help me track down the cause of the apparent undershoot in the pre-mini baseline in the average mini.
# Fixed undershoot by aligning on mid-rise, but still having trouble excluding some poor minis. 
# Will add a linear fit to whole mini to compare with 3-part fit.
#
#Current version: mini4Mar29
#
# Found a bug by which minis are double counted if there are two points detected within a brief period. 
#
# will exclude in measure mini
#
# Also, adding codes to indicate why minis are being rejected, should aid adjusting parameters if desired. Putting a threshold on integral, after
# examining a lot of data it seems integral <150 are almost all bad and <250 are at least half bad. So set at 250 for now. 
#
# Adding separate plotting routine to save the plotly plots of individual trials. The purpose
# is to color the good and bad minis within the analyzed wave for easier trouble shooting of what is being
# accepted and rejected. Also adding an option to plot waves as negative going (actual polarity) vs. always showing them positive
#
# This is mini4Feb11.R 
#
# Minor upgrade to make folder naming expectations more flexible
# Latest Version Nov 7, now called mini4Nov4.R
# fixed exclusion criteria
# fixed bug which crashed if last mini was cut off
# 7/22/22 mini4.r 
# modified from mini3, to batch analyze all cells.
# Assumes directory structure with each cell in its own directory within the user chosen directory.
# Bad traces are detected from errors and put in a subdirectory.
# Plots of traces with detected minis are saved to HTML savewidget (requires pandoc)
# removed single file mode and flags to permit getting or not getting data and detecting or not detecting minis
#
# TO DO:
# 1. need to assess goodness of fit and reject minis as needed, probably can use broom to
#    save tidy output from fits and then measure significance and residuals. 
# 2. savewidget saves single file, but also saves many other files, not sure if it is self contained, need to test
# 3. ideally would have shiny UI for setting parameters--this could be separate program that then writes parameter 
#       to a file.
# 4. summary histograms--both for individual cells and sampled across cells, the latter will also need
#       an interface for assigning conditions to cells--this could be added to the shiny UI.
#####################################################################

library(data.table)
library(tidyverse)
library(RcppRoll)
#library(Rfast)
library(IgorR)
library(plotly)
library(htmlwidgets)
library(DT)
library(minpack.lm)
library(gsignal)
library(RColorBrewer)
library(broom)
if (!(is.null(options("viewer")))) {
  backupOptions <- options()
  options(viewer=NULL)                  # this is important for viewing plotly plots in separate windows
}

Sys.setenv(RTUDIO_PANDOC="~.local/share/pandoc")
#************************************************************************
#*
#* Functions NOW moved to genUtils.R and ephysUtils.R; these must be loaded and sourced.
#* Could include the source statement here, but more flexible to not insist on a specific location
#*
#************************************************************************


#************************************************************************
#*
#* measure minis found by template match
#*
#************************************************************************

measurePks <- function(titl, tDat, riseLocs, dPar) {
  #
  # locations of minis detected by template match are in riseLocs, but are offset by
  # baseline period in template, and by any difference in kinetics from template.
  # 1) find peak within max expected period
  # 2) fit line to baseline and to rise and define start from intersection
  # 3) fit decay
  # 4) store mini for later calculation of average mini.
  #

  minRise = dPar$minRise # start has to be at least 1 pt before late rise, or reject
  minNextRise = dPar$minNextRise
  blSubt = dPar$blSubt
  midRisePt <- dPar$templBase + floor(dPar$maxTTPS/2) # align mid rises of extracted minis at this point (e.g. 50 +40/2 = 70)
  
  #
  # setup data table of results returned each datafile
  #
  npks <- length(riseLocs)
  pkLoc <- vector(mode="integer", length=npks)
  TTP <- vector(mode="integer", length=npks)
  pkAmp <- vector(mode="double", length=npks)
  rise20 <- vector(mode="integer", length=npks) 
  rise80 <- vector(mode="integer", length=npks)
  rise10 <- vector(mode="integer", length=npks)         # not sure whether we will want to use 20-80 or 10-90, so compute both
  rise90 <- vector(mode="integer", length=npks)
  mStrts <- vector(mode="integer", length=npks)
  intgrl <- vector(mode="double", length=npks)
  baseIncpt <- vector(mode="double", length=npks)
  baseIp <- vector(mode="double", length=npks)
  baseSlope <- vector(mode="double", length=npks)
  baseSp <- vector(mode="double", length=npks)
  baseResid <- vector(mode="double", length=npks)
  baseMn <- vector(mode="double", length=npks)
  riseIncpt <- vector(mode="double", length=npks)
  riseIp <- vector(mode="double", length=npks)
  riseSlope <- vector(mode="double", length=npks)
  riseSp <- vector(mode="double", length=npks)
  riseResid <- vector(mode="double", length=npks)
  riseConf1 <- vector(mode="double", length=npks)
  riseConf2 <- vector(mode="double", length=npks)
  
  dkAmp <- vector(mode="double", length=npks)
  dkAmpP <- vector(mode="double", length=npks)
  dkTau <- vector(mode="double", length=npks)
  dkTauP <- vector(mode="double", length=npks)
  dkResid <- vector(mode="double", length=npks)
  totResid <- vector(mode="double", length=npks)
  totFit <- vector(mode="double", length=npks)
  linRatio <- vector(mode="double", length=npks)
  linIncpt <- vector(mode="double", length=npks)
  linIp <- vector(mode="double", length=npks)
  linSlope <- vector(mode="double", length=npks)
  linSp <- vector(mode="double", length=npks)
  linResid <- vector(mode="double", length=npks)
  pAmp <- vector(mode="double", length=npks)
  dAmp <- vector(mode="double", length=npks)
  pRise <- vector(mode="double", length=npks)
  dRise <- vector(mode="double", length=npks)
  pDK <- vector(mode="double", length=npks)
  dDK <- vector(mode="double", length=npks)
  
  xclude <- vector(mode="integer", length=npks)    # should be initialized to 0, so set to + to exclude
  
  sweepRslts <- list()                # item 1 is data.table of results, item 2 is array of extracted minis
  bl <- dPar$templBase                # baseline of template mini (# of samples)
  minBl <- bl/2                       # rejecting too many minis because the starts are <50
  maxTTP <- bl + dPar$maxTTPS         # latest allowable peak including baseline for good minis
                                      # note, this is different from maxPkFound, the latest we look for a peak
  dkLen <-dPar$maxDkS
  mLen <- maxTTP + dkLen              # length of a mini
  dLen <- nrow(tDat)                  # length of the data
  bVerbose <- dPar$bVerbose
  
#
# extracted minis have at least templBase points of baseline followed by the full mini (defined by where rise intersects
# baseline), followed by additional points to pad to the maximum allowable length mLen. The minis are aligned so the midrise
# falls on midRisePt
#
  xMin <- matrix(ncol=npks, nrow=mLen)
  xMinf <- matrix(ncol=npks, nrow=mLen)
  xFit <- matrix(ncol=npks, nrow=mLen)
  
  pad = trunc(dPar$mPkPts/2)                # padding needed for rolling averages

  endSamp <- 0                              # end of portion of datafile displayed at one time, starts at 0
print(paste("#peaks=",length(riseLocs)))

  for (ipk in 1:npks) {
    print(paste0("ipk=",ipk))
#
# first find peak. Use rolling average to smooth, but this requires padding (e.g. 3 pt rolling avg deletes
#  1 pt on either side), so add back the pad to retrieve the correct index in the original vector.
#
    strt <- riseLocs[ipk]
    if (ipk > 1) {
      if (strt - riseLocs[(ipk-1)] <= minNextRise) { # skip it if it overlaps rise of last peak
        print(paste0("overlaps"))
        xclude[ipk] <- 1              # exclude code 1 = overlaps
      }  
    }
    arise <- roll_mean(tDat$yf[strt:(strt+maxTTP)],dPar$mPkPts) # this should contain baseline + rise and peak
    baseAmp <- mean(arise[1:(minBl)])               # take mean of first minBl points as initial baseline and baseline subtract
    arise <- arise - baseAmp
    pkLocInRise <- which.max(arise[1:maxTTP])
    pkLoc[ipk] <- pkLocInRise + strt + pad - 1              # since strt and pkLocInRise are both 1-based, we are double counting, so subtract 1.
    pkA <- arise[pkLocInRise]
#browser()
#
# copy y values into a vector, find mid rise and then extract mini, aligned on mid rise
#
print(paste("strt=", strt,"pkA1=",pkA,"pkLocInRise=",pkLocInRise, "pkLoc=", pkLoc[ipk]))
    midAmp <- 0.5 * pkA
    midRise <- detect_index(.x=arise[1:pkLocInRise], .f=function(x) x <= midAmp, .dir="backward")       # look backward to find 50% of pk
#
# error in midrise with sharp minis so reverting to single look back from pk
#
    mrInSw <- midRise + strt + pad - 1 # location of mid rise in WHOLE SWEEP
    print(paste("midRise=",midRise,"mrInSw=",mrInSw, "midRisePt=",midRisePt, "mLen=",mLen))
    xBeg <- mrInSw-midRisePt+1
    xNd <- mrInSw + (mLen - midRisePt)
    if ((xBeg > 0) & (xNd < dLen)) {                   # if there is room to extract the mini without overrunning
print(paste("extracting mini from =",xBeg, ":", xNd))
      xMin[1:mLen, ipk] <- tDat$yf[xBeg:xNd]    # storing median filtered wave
      xMinf[1:mLen, ipk] <- tDat$yff[xBeg:xNd]  # storing low pass filtered wave
      baseAmp2 <- mean(xMin[1:bl, ipk])   # correct baseAmp with full baseline period
      xMin[1:mLen, ipk] <- xMin[1:mLen, ipk] - baseAmp2
      xMinf[1:mLen, ipk] <- xMinf[1:mLen, ipk] - baseAmp2
      vMin <- xMin[, ipk]    # single vector versions
      vMinf <- xMinf[, ipk]
      baseMn[ipk] <- round((baseAmp+baseAmp2), 1)
#
#   refind the peak after changing alignment, baseline subtracting
#
#      
      rollMini <- roll_mean(vMin,dPar$mPkPts)
      pkLocInMin <- which.max(rollMini[1:(dPar$maxPkFound)])
      TTP[ipk] <- pkLocInMin
#print(paste("pkLocInMin=",pkLocInMin, "maxTTP=", maxTTP))
      if (bVerbose) print(paste0("TTP=",TTP[ipk]))
#if (ipk==2) browser()
      pkA <- rollMini[pkLocInMin]
      pkAmp[ipk] <- round(pkA,1)
      print(paste("pkLocInMin=",pkLocInMin, "pad=",pad, "rollMini[pkLocInMin]=",rollMini[pkLocInMin], "pkA=",pkA))
    } else {
      xclude[ipk] <- 12                   # not enough room at beginning or end of data.
    }
print(paste("midrise in arise=",midRise, "midAmp=", midAmp,"xmin=",xMin[midRisePt,ipk]))
#if(ipk==5) browser()
# 
# now find the other rise levels
#
#
# NOTE: not error checking to make sure >80% and <20% values exist. A strategy based on match(TRUE, revRise > amp) returns NA if 
# no true elements occur, but which.max returns the first element if no true elements occur.
#
# ALSO NOTE: fitting rise in filtered wave!!!
#
    if (xclude[ipk] <=0) {
      lateRise8 <- detect_index(.x=vMinf[1:pkLocInMin], .f= function(x) x < 0.8*pkA, .dir="backward")       # look backward to find 80% of pk
      lateRise9 <- detect_index(.x=vMinf[1:pkLocInMin], .f= function(x) x < 0.9*pkA, .dir="backward")       # look backward to find 90% of pk
      mStrt <- 0 # default for printing out if we don't find an actual start
      if ((lateRise8 == 0) | (lateRise9 == 0)) {
        xclude[ipk] <- 2                                      # exclude code 2 = can't find lateRise
        earlyRise1=0
        earlyRise2 = 0
      } else {
        earlyRise1 <- detect_index(.x=vMinf[1:lateRise8], .f= function(x) x < 0.1*pkA, .dir="backward")      # same for 10% of pk
        earlyRise2 <- detect_index(.x=vMinf[1:lateRise8], .f= function(x) x < 0.2*pkA, .dir="backward")      # same for 20% of pk
      }
      if ((earlyRise1 == 0) | (earlyRise2==0)) {
        xclude[ipk] <- 3                                      # exclude code 3 = can't find earlyRise
      } else {
        mStrt <- detect_index(.x=vMinf[1:earlyRise1], .f= function(x) x < 0.05*pkA, .dir="backward")        # start of mini defined as 5% of peak
        if (mStrt==0) xclude[ipk] <- 4                        # exclude code 4 = can't find mStrt
      }
      if (bVerbose) print(paste("rises=",earlyRise1, earlyRise2, lateRise8, lateRise9, " mStrt=",mStrt, "pkA=",pkA))
    } # end if xcluding
    if (xclude[ipk] <=0)  {                                 # if not excluding
      rise20[ipk] <- max(0, (earlyRise2 - mStrt))            # subtracting starting index from times
      rise80[ipk] <- max(0, (lateRise8 - mStrt))
      rise10[ipk] <- max(0, (earlyRise1 - mStrt))            # subtracting starting index from times
      rise90[ipk] <- max(0, (lateRise9 - mStrt))
      if ((mStrt > minBl) & (rise80[ipk] > minRise)) {      # have to have at least a couple of points of rise to fit
        mStrts[ipk] <- mStrt                                # NOTE: mStrts can be < bl
      } else {
        xclude[ipk] <- 5                                    # exclude code 5 = mStrt not > minBL or rise80 not > minRise
        if (bVerbose) print(paste("rejecting since mStrts[]ipk", mStrts[ipk], " is < minBl", minBl))
      }
    } else {
      if (bVerbose) print(paste("rejecting mini code=",xclude[ipk]))
    }
#
# For now, rejecting minis that don't have mStrt > bl or if early or late rise not found
# NOTE: times are offset by strt which is the same as riseLocs[ipk]
#
# now fit decay, get params--offset, so all values are positive, to avoid log problems
#

#
# fit baseline and rise, since mStrt now known
# The way this works is lm returns a fit object and the broom functions augment and tidy
# return dataframes with fitted values, residuals etc. (augment) and fit parameters (tidy)
#
# HERE FITTING ORIGINAL WAVE
#
# NOTE: baseline includes points up to mStrt-1 (mStrt not included)
#       rise includes mStrt and thePk
#       decay is from thePk to thePk + dkLen-1; in other words thePk is both last point of
#         rise and first point of decay. But have to not double save this point
#
# First do linear fit to whole mini to compare to 3-part fit
# For sharp minis, changing this to exclude baseline
#
    if (xclude[ipk] <=0) {
    #
      theBase <- data.table(x = 1:bl, y=vMin[1:bl])
      bFit <- lm(y ~ x, theBase)
      bAug <- augment(bFit, theBase)
      bTid <- tidy(bFit)
      baseIncpt[ipk] <- round(bTid$estimate[1],1)
      baseIp[ipk] <- signif(bTid$p.value[1],2)
      baseSlope[ipk] <- round(bTid$estimate[2],2)
      baseSp[ipk] <- signif(bTid$p.value[2],2)
      baseResid[ipk] <- round(sum(abs(bAug$.resid)),1)
# now rise
      x1 <- earlyRise1; x2 <- lateRise9                             # CHANGING to fitting 10-90 rise, not 5 to peak
      theRise <- data.table(x=x1:x2, y=vMin[x1:x2])
      riseDur = x2 - x1 + 1
      if (all(is.na(theRise$y)))  {
        xclude[ipk] <- 6                                # exclude code 6 = can't fit line to rise
        if (bVerbose) print(paste("Unable to fit rise from ",x1," to ", x2))
      } else {
        riseFit <- lm(y ~ x, theRise)
        rAug <- augment(riseFit, theRise)
        riseResid[ipk] <- round(sum(abs(rAug$.resid)),2)
        rTid <- tidy(riseFit)
        riseIncpt[ipk] <- round(rTid$estimate[1],1)
        riseIp[ipk] <- signif(rTid$p.value[1],2)
        riseSlope[ipk] <- round(rTid$estimate[2],2)
        riseSp[ipk] <- signif(rTid$p.value[2],2)
      }
#if (ipk==5) browser()   # & (titl=="C10T2")
#print(paste("riseIncpt=",riseIncpt[ipk],"riseSlope=",riseSlope[ipk],))
#
# Since length of mini is maxTTP (inc baseline) + max decay, and we are fitting max decay, we will overrun if
# TTP is > maxTTP. Could fit a subset if needed, but might as well just not fit, since mini is "bad" anyway.
# so don't proceed e.g. with augment if fit is bad or TTP > maxTTP
#
      mDk <- data.table(t=1:dkLen, y=vMinf[pkLocInMin:(pkLocInMin+dkLen-1)])       # 
      dFitTst <- try(nlsLM(y ~ expDecay1(t,amp,tau), data=mDk, start=list(amp=pkAmp[ipk], tau=dPar$tauDkT)))
      if((inherits(dFitTst, "try-error")) | (TTP[ipk] > maxTTP)) {
        if (bVerbose) print(paste("error fitting decay"))
      } else {  
        dAug <- augment(dFitTst, mDk)
        dkResid[ipk] <- round(sum(abs(dAug$.resid)),2)
        dTid <- tidy(dFitTst)
        dkAmp[ipk] <- round(dTid$estimate[1],1)
        dkAmpP[ipk] <- signif(dTid$p.value[1],2)
        dkTau[ipk] <- round(dTid$estimate[2],2)
        dkTauP[ipk] <- signif(dTid$p.value[2],2)
        if (bVerbose) print(paste("dkAmp=", dkAmp[ipk], "dkAmpP=", dkAmpP[ipk], "dkTau=", dkTau[ipk]))

        if((is.nan(sum(bAug$.fitted))) | (is.nan(sum(rAug$.fitted))) |
            (is.nan(sum(dAug$.fitted)))){
          print(paste("mStrt=",mStrt,"pkLocInMin=", pkLocInMin, "riseDur=",riseDur,"dkLen=",dkLen,"mLen=",mLen))
          browser()
        }
        xFit[1:mLen, ipk] <- xMin[1:mLen, ipk]        # copy mini into fit, a few points will remain from original data, rest will be replaced by fits
        xFit[1:bl, ipk] <- bAug$.fitted[1:bl]
        xFit[x1:x2, ipk] <- rAug$.fitted[1:riseDur]
        xFit[pkLocInMin:(pkLocInMin+dkLen-1), ipk] <- dAug$.fitted[1:(dkLen)]
        totFit[ipk] <- round(sum(xFit, na.rm = T),1)
        totResid[ipk] <- round(sum(baseResid[ipk],riseResid[ipk],dkResid[ipk],na.rm = T), 1)
        riseDur <- pkLocInMin - bl                     # duration of rise is from end of baseline to peak
        minBeg <- pkLocInMin - riseDur + 1             # this is the first point of the rise (i.e. after baseline)
# now linear fit to rise and initial decay to compare goodness of fit
        tau4 <- round(4*dkTau[ipk])
        fitBeg <- max(1,(minBeg-tau4))
        fitNd <- min((pkLocInMin + tau4), length(vMin)) # end of portion of decay we are fitting--don't overrun the end.
        if ((fitBeg >0) & (fitNd > fitBeg)) {             # make sure there are non-zero vals
          riseDK <- data.table(x = 1:(fitNd-fitBeg+1), y= vMin[fitBeg:fitNd]) # use Y values from start of mini to 4 * decay tau, or as much as will fit
          lFit <- lm(y ~ x, riseDK)
          lTid <- tidy(lFit)
          lAug <- augment(lFit, riseDK)
          linIncpt[ipk] <- round(lTid$estimate[1],1)
          linIp[ipk] <- signif(lTid$p.value[1],2)
          linSlope[ipk] <- round(lTid$estimate[2],2)
          linSp[ipk] <- signif(lTid$p.value[2],2)
          linResid[ipk] <- round(sum(abs(lAug$.resid)),1)
        print(paste("linResid=",linResid[ipk]))
#        linRatio[ipk] <- round((linResid[ipk]/totResid[ipk]),2)
          intgrl[ipk] <- round(sum(vMin[fitBeg:fitNd]))            # compare residual of linear fit to integral over same range    
          if (intgrl[ipk] > 0) linRatio[ipk] <- round((linResid[ipk]/intgrl[ipk]),2)
        } else {
          print(paste("can't do linear fit"))
        }
        pdlst <- miniStat(vMinf, minBeg, pkLocInMin, dkTau[ipk], 3)
        if (pdlst$err==0) {
          pd <- pdlst$pd
          pAmp[ipk] <- pd$pWelchT[1]
          dAmp[ipk] <- pd$cohensD[1]
          pRise[ipk] <- pd$pWelchT[2]
          dRise[ipk] <- pd$cohensD[2]
          pDK[ipk] <- pd$pWelchT[3]
          dDK[ipk] <- pd$cohensD[3]
        } else {
          print(paste("error from miniStat--too long a decay tau?"))
        }
#browser()
        if (dPar$seeEach) {
          theMini <- data.table(x=1:mLen, y=xMin[1:mLen,ipk], yfit=xFit[1:mLen, ipk])
          p <- plot_ly(theMini, x = ~x, y = ~y, type="scatter", mode="markers")
          p <- p %>% add_trace(y = ~yfit, mode="lines")
          print(p)
          readline(prompt = '')
        }  
      } # if no error fitting decay
    } # end if not rejecting
  } # end for minis in trace
#
# Now exclude minis that do not meet criteria
# EXLUSION CRITERIA (copied from beginning of Main above:
# maxBslope = 0.1  # max  of abs(baseline slope) 
# minPk = 5      
# maxTTPS <- 30    # max allowable mini rise duration in samples: (=3 ms at 10 kHz)
# minRiseSlope <- 1
#
# exclude code 1 = overlaps
# exclude code 2 = can't find lateRise
# exclude code 3 = can't find earlyRise
# exclude code 4 = can't find mStrt
# exclude code 5 = mStrt not > minBL or rise80 not > minRise
# exclude code 6 = can't fit line to rise
# exclude code 11 = peak occurs too late
# exclude code 12 = not enough room at beginning or end of data.

  xclude[which((abs(baseSlope) > dPar$maxBslope)&(xclude==0))] <- 8   # exclude code 8 = baseline slope too high
  xclude[which((pkAmp < dPar$minPk)&(xclude==0))] <- 9                # exclude code 9 = amplitude too small
  xclude[which((riseSlope < dPar$minRiseSlope)&(xclude==0))] <- 10    # exclude code 10 = rise slope too small
  xclude[which((TTP > maxTTP)&(xclude==0))] <- 11                     # exclude code 11 = time to peak too late
  xclude[which((riseSp > dPar$maxRiseSp)&(xclude==0))] <- 13          # exclude code 13 = poor fit to rise
#  xclude[which((linRatio < dPar$minLinRatio)&(xclude==0))] <- 14     # exclude code 14 = poor overall fit, relative to straight line
  
#  xclude[which((dDK < dPar$minDDK) & (xclude==0))] <- 14              # exclude code 14 = poor effect size for decay slope
# not excluding for decay and rise for now
#  xclude[which((dRise < dPar$minDrise) & (xclude==0))] <- 15          # exclude code 15 = poor effect size for rise slope
  xclude[which((intgrl < dPar$minIntgrl) & (xclude==0))] <- 7         # exclude code 7 = low integral
  #browser()
  sweepRslts$res1 <- data.table(riseLocs, pkLoc, TTP, pkAmp, rise10, rise20, rise80, rise90, mStrts, intgrl, baseIncpt,
      baseIp, baseSlope, baseSp, baseResid, baseMn, riseIncpt, riseIp, riseSlope, riseSp, riseResid, dkAmp, dkAmpP, dkTau, dkTauP, dkResid, totResid, totFit,
      linRatio, linIncpt, linIp, linSlope, linSp, linResid, pAmp, dAmp, pRise, dRise, pDK, dDK, xclude)
  sweepRslts$xMin1 <- xMin
  sweepRslts$xMin1f <- xMinf
  sweepRslts$xFit1 <- xFit
  return(sweepRslts)
}

################################################################
#************************************************************************
#*
#* Brute force template matching mini detection function
#*
#************************************************************************

detectPks <- function(tDat, templ, dPar) {
  
  #
  # for now match at each sample, but could do only some (e.g. every other) since match is likely to be
  # pretty good even if off by 1 sample.
  #
  # To start just focus on detection, can use extraction and fitting code from detectPks above, once detection
  # works
  #
  tttp <- dPar$ttpTempl
  mLen <- dPar$templBase + dPar$maxTTPS + dPar$maxDkS
  tLen <- length(templ)     # template length
  dLen <- nrow(tDat)        # data length
# make indices
  nTempl <- dLen-tLen+1  # #of templates that fit into data
  indV <- vector(mode="integer", length=(tLen*nTempl))
  for (i in 0:(nTempl-1)) {
    indV[(1+tLen*i):(tLen*(i+1))] <- (1:tLen)+i
  }
#browser()
  
#
# now template match
#
#Put data in matrix; each template length of data is in a column, each row is shifted by 1
# Then normalize each column to the point corresponding to the peak of the template
#
#  browser()  
  
  dMat<- matrix(data=tDat$yff[indV],nrow=tLen, ncol=nTempl)
#
# now linearly regress template against matrix, yields a vector of coefficients (amp and intcpt)
#
  fm <- lm(dMat ~ templ)
  coefDT <- data.table(t(coef(fm)))     # need to transpose since rows are (Intercept) and templ
  setnames(coefDT, '(Intercept)', 'intcpt')
#
# now store threshold crossings and return
# We have runs of values above threshold, for each run, we need to find the best one
# To detect runs, subtract a lagged version and then use rle
#
  minVal = dPar$SDthresh * sd(tDat$yf)
  aboveTh <- which(coefDT$templ>=minVal)
  shiftDiffTh <- aboveTh - data.table::shift(aboveTh, type="lag", fill=(aboveTh[1]-1))
  shRLE <- rle(shiftDiffTh)
#  
# for each run, get starting index into aboveTh, which in turn is an index into the trace
# 
  runStrt <- cumsum(shRLE$lengths)
  bRuns <- which(shRLE$values==1)      # all runs have values of 1 
  run1End <- runStrt[bRuns]            # runs end on the corresponding cumulative index
  # the start is the length of the run earlier + 1
  run1Strt <- run1End - shRLE$lengths[bRuns] + 1
  npks <- length(run1Strt)
  if (npks >0) {
#
# runs end on each element of run1End and start on run1Strt
#
    riseLocs <- vector(mode="double", length = npks)
#browser()
    for (ir in 1:npks) {
      maxFitInRun <- which.max(coefDT$templ[aboveTh[run1Strt[ir]:run1End[ir]]]) # should find the max
      riseLocs[ir] <- aboveTh[run1Strt[ir]] + maxFitInRun - 1 
      if ((riseLocs[ir] + mLen) >= dLen) riseLocs[ir] <- 0      # reject if not enough room for a mini
    }
  } else {
      riseLocs = 0
  }
  riseLocs <- riseLocs[riseLocs>0]
#
# riseLocs should actually be the start location at this point (or close to it) since this is where the
# the fit to the template has the maximum amplitude and the template also includes baseline
# riseLocs is in samples, since no time vector included
  return(riseLocs)
}
#************************************************************************
#*
#* plotTrl
#*
#************************************************************************
plotTrl <- function(tDat, titl, nMini, dPar, mRes, outdir) {
 
  pkLocs <- mRes$pkLoc 
  bl <- dPar$templBase                # baseline of template mini
  maxTTP <- bl + dPar$maxTTPS         # latest allowable peak including baseline, each mini starts at this point before peakloc
  mLen <- maxTTP + dPar$maxDkS        # length of each mini
  tDat$"colr" <- 1        # set points to 1 of 3 colours: black=1 (default), red=2 rejected minis, green=3: accepted minis
  nMini <- length(pkLocs)
  xcl <- mRes$xclude
  nRejected <- sum(xcl > 0)       # number xcluded
  nAccepted <- nMini-nRejected
  rPks <- vector(mode="integer", length=(nRejected))
  aPks <- vector(mode="integer", length=(nAccepted))
  rPks <- pkLocs[which(xcl > 0)]  
  aPks <- pkLocs[which(xcl <= 0)]
  txtVals <- as.character(xcl[which(xcl>0)])      # take the non-zero exclusion codes and turn them into characters
  
#  browser()
  fig <- plot_ly(tDat, x=~t, y=~yff, color=I("black"), name='filt', type='scatter', mode='lines')
  fig <- fig %>% add_lines(y=~yf, name='orig ', linetype='dot', alpha=0.2)
  fig <- fig %>% add_trace(x=~aPks/dPar$sr, y=~yf[aPks], color=I("green"), name='aPks', mode='markers')
  fig <- fig %>% add_trace(x=~rPks/dPar$sr, y=~yf[rPks], color=I("red"), name='aPks', mode='markers', hoverinfo="text", hovertext=txtVals)
#  print(fig)
  l <- layout(fig, title = titl)
#print(paste0("outdir=",outdir," titl=",titl))
  outf <- paste0(outdir,titl,"index.html")
#
# htmlwidgets fails if not connected to the internet
  htmlwidgets::saveWidget(as_widget(partial_bundle(l)), outf, selfcontained = T)
#
# saveWidget is supposed to save a self contained file without also saving a folder of other stuff.
# the index file is self contained, but the extra folder is still created, so just get rid of it.
#
  outfold <- str_replace(outf,"index.html","index_files")
  unlink(outfold, recursive = T)

}
#************************************************************************
#*
#* Main file loop
#*
#************************************************************************

minis <- function(defSR=10000, defSWdur=10) {

  sweepRslts <- list()
  bVerbose = T

# filtering parameters
  mFiltPts = 3   # keep this odd so we can always take the center as peak or trough
  cutOff = 1000     # cutoff frequency for lowpass butterworth filter
  superLow = 10 # cutoff for filtered version used to remove baseline fluctuation
  bMedFilt = T   # flag for turning on/off median filtering

# detection paramters
# time parameters are optimized for sampling rate of 10Kz, so adjust by multiplying by sr/10000 in case sr is altered
  sr = defSR     # correct SR after reading file if need be
  swDur = defSWdur     # duration of data in file, correct if needed
  minPk = 5      #
  aScale = 1e12  # convert amps to pA
  mPol = -1      # detect negative going minis (scale converts to positive going)
  mPkPts = 3     # detect peak in original, so average over 3 pts
  fitSealDur = 0.1  # fit only initial 100 ms after peak of seal test
#################################################################################################
#------------------------------------------------------------------------------------------------
  
  sealPad = 0.03    # period for relaxation of return to baseline after seal test
  SDthresh = 2  # use a threshold of two SD of median filtered data.                 TESTING 1.5
  maxBslope = 0.15  # max  of abs(baseline slope) 
  numTauDK = 1    # fit one exp to decay (other allowable value is two)
  minRise = 1    # have to be at least 2 points to fit rise
  maxRiseSp = 0.02 # probability of rise fit must be less than this
  minNextRise = 30 # exclude overlapping rises within 30 samples
#  minIntgrl = 250 # empirically determined looking at data
  minIntgrl = 200 # empirically determined looking at data
  minDDK = 0.5    # min effect size for decay slope < than baseline
  minDrise = 0.5 # min effect size for rise slope > than baseline
  detrendMode = "linear" # "filter" to subtract a very low freq high pass version, or "linear" to just do a linear detrending
#
# template and mini extraction parameters. Strategy is to just calculate fit for period of template.
# template will include moderately fast rise and decay and no baseline. Extracted mini and template must be aligned on peak.
#
#
  templBase <- 50   # baseline of 5 ms at 10KHz: used for template and minimal baseline of extracted minis
  blSubt = templBase  # baseline subtract mean of this baseline period. Use blSubt = 0 to disable baseline subtraction.
  maxTTPS <- 50    # max allowable mini rise duration in samples: (=5 ms at 10 kHz) for good mini
  minRiseSlope <- 0.12 # slope of rise at least 0.12 pA/sample = 1.2pA/ms, so will rise to almost minPk amp threshold of 5 pA in maxTTPs (max rise duration) of 5 ms
  minLinRatio <- 2.0  # residuals from simple linear fit at least this many times summed residuals from 3-part fit 
#  minLinRatio <- 1.8  # residuals from simple linear fit at least this many times summed residuals from 3-part fit 
  maxTDk <- 80     # max tau-decay in samples  (10 ms at 10kHz)
  maxPkFound <- templBase + maxTTPS + maxTDk    # latest we will look for peak including baseline--one max decay tau past max peak
  
#-------------------------------------------------------------------------------------------------
##################################################################################################
    
  nTauDk <- 3       # extract candidate mini with enough room for this # * maxTDk
  maxDkS <- nTauDk*maxTDk   # total maximum decay duration for extracted minis
  mLen <- templBase + maxTTPS + maxDkS      # length of mini is max baseline + ttp + max decay from peak. 
#
# display parameters
#
  seeEach = F
  bDoRinVhold = T
  dspSamps = 2 * sr   # duration of displayed trace in samples (individual mini display durations are fixed); 0 means don't display
#
# make mini template
#
  ttpTempl <- 10     # ttp for template in samples (assume 10khz): 1 ms * sr; NOTE last pt is pk, so rise occurs over ttp-1 increments
  tauDkT <- 40      # tau of 4 ms * sr for template
#  tauDkT <- 15      # tau of 1.5 ms for template.            TRYING MUCH SHARPER MINI FOR DEREKS DATA
  nTauTempl <- 5    # make template with decay of 5* tauDkT
  templDKlen <- tauDkT*nTauTempl
  tlen <- templBase + ttpTempl + templDKlen      # total duration of template=rise (inc peak) + decay
  templ <- vector(mode="numeric", length=tlen)
  templ[1:templBase] <- 0
  templ[(templBase+1):(templBase+ ttpTempl)] <- seq(from=0, to=1, by = (1/(ttpTempl-1)))  # assumes a peak amp of 1
  t <- seq(1, (templDKlen), by=1)                  # time in samples, starting after rise
  templ[(templBase+ttpTempl+1):tlen] <- exp(-t[1:templDKlen]/tauDkT)
#
#
# put needed parameters into list
#

  dPar <- list(
    "seeEach"     = seeEach,         # whether to display each mini and fits
    "sr"          = sr,              # sampling rate in Hz
    "swDur"       = swDur,           # duration (s) of one sweep of data (i.e. 1 file)
    "minPk"       = minPk,
    "minRise"     = minRise,
    "minNextRise" = minNextRise,
    "minIntgrl"   = minIntgrl,
    "minDDK"      = minDDK,
    "minDrise"    = minDrise,
    "minLinRatio" = minLinRatio,     # ratio of residuals for straight line fit vs. 3-fold baseline, rise, decay fit.
    "cutoff"      = cutOff,          # for butterworth filter
    "mPkPts"      = mPkPts,          # number of samples to average mini peak over
    "maxTTPS"     = maxTTPS,         # max allowable time to peak (not inc baseline) for good mini
    "maxPkFound"  = maxPkFound,      # latest we will look for a peak
    "maxTDk"      = maxTDk,          # max allowable tau decay
    "maxDkS"      = maxDkS,          # max decay phase in samples of extracted minis
    "tauDkT"      = tauDkT,          # template decay constant
    "SDthresh"    = SDthresh,        # threshold for template fit to mini based on SD of median filtered signal
    "templBase"   = templBase,       # length of baseline for extracted minis and template
    "blSubt"      = blSubt,          # baseline subtraction window (or 0 to disable subtraction)
    "maxBslope"   = maxBslope,       # max abs value of baseline slope
    "ttpTempl"    = ttpTempl,        # TTP used for template
    "minRiseSlope"= minRiseSlope,    # minimum allowable rise slope
    "mWaveStrt"   = 0L,              # sample number at which mini portion of wave starts, calculated below
    "dlen"        = 0L,              # length in samples of mini containing portion of data, calculated below 
    "plotInv"     = F,               # plot data in original polarity for trial-long plots
    "detrendMode" = detrendMode,     # "linear fit" or "filter"
    "bVerbose"    = bVerbose)
#
# Get data
#
  inFold <- rstudioapi::selectDirectory(caption = "Select Condition Folder", label = "Select")
  prefix <- "Cell"; nCells <- 0
  

  if (is.null(inFold)) break
#browser()  
  dirList <- list.dirs(path=inFold, recursive = F)
  if (sum(str_detect(dirList,"results_")) > 0) { 
    reply <- getUserInput(paste("This looks like a single cell folder, y or Y to proceed in single cell mode, or anything else to pick again"))
    if ((reply=="Y") | (reply=="y")) {
      nCells == 1
    } else {
      minis()
    }
  } # if looks like single cell folder
  if (sum(str_detect(dirList,regex(prefix,ignore_case=T))) == 0) {
    while (nCells==0) {
      nCells <- sum(str_detect(dirList,regex(prefix,ignore_case=T)))
      if (nCells==0) {
        prefix <- getUserInput(paste("no enclosed folders with ", prefix, " enter different cell folder string (case ignored):  "))        # default varType is string, so not needed
        nCells <- sum(str_detect(dirList,regex(prefix,ignore_case=T)))
      }
    }
  } else {
    nCells <- sum(str_detect(dirList,regex(prefix,ignore_case=T)))
  }
  cellFolds <- dirList[which(str_detect(dirList,regex(prefix,ignore_case=T)))]
  conditionFold <- inFold
#browser()
#
# test first file name for 4 digit format for cells, trials
#
# data file parameters
  sep <- "_"
  nCellDig <- 4
  nTrialDig <- 4
  extPatt <- ".ibw$"     # regex for .ibw at end of string
  cellStr <- vector(mode = "character", length = nCells)
  cellGdMins <- vector(mode = "integer", length = nCells)
  cellAmp <- vector(mode = "numeric", length = nCells)
  
  
  for (iCell in 1:nCells) {
    inFold <- cellFolds[iCell]                            # using inFold, so I don't have rewrite code below
    cellStr[iCell] <- basename(cellFolds[iCell])
    flist <- list.files(path=inFold, pattern=extPatt)     # stored prior inFold in conditionFold
    fName <- flist[1]
#
#   Create folder for results
#
    outdir1 <- paste0(inFold,"/results_",Sys.Date(),"/")
    bExit=F
    resNum = 0
    while (bExit==F) {
      if (dir.exists(outdir1)==FALSE) {
        outdir <- outdir1
        dir.create(outdir)
        bExit=T
      } else {
        resNum = resNum+1
        outdir1 = paste0(inFold,"/results_",Sys.Date(),letters[resNum],"/")
        if (resNum >25) {
         stop("too many results files")
        }
      }
    } # end while not exiting
#
# write detection parameters
#
    
    nFiles <- length(flist)
    if (bVerbose) print(paste(nFiles," data files"))
    trial <- vector(mode = "integer", length = nFiles)
    titl <- vector(mode = "character", length = nFiles)
    rins <- vector(mode = "numeric", length = nFiles)
    vHold <- vector(mode = "numeric", length = nFiles)
    nMini <- vector(mode = "integer", length = nFiles)
#    outRes <- vector(mode = "character", length = nFiles)      # these were paths for results files, don't think using now
#    outMin <- vector(mode = "character", length = nFiles)
#    outMinF <- vector(mode = "character", length = nFiles)
#    outFit <- vector(mode = "character", length = nFiles)
    rms <- vector(mode = "numeric", length = nFiles)
    ndPre <- nchar(prefix) + nchar(sep)
    ndCell <- ndPre + nCellDig
    for (iFil in 1:length(flist)) {
      fName <- flist[iFil]
#
# test fName and correct defaults if needed
#
      if(str_detect(fName,sep)==F) {
        sep <- getUserInput(paste("change sep from ", sep," to: "))
      }
      if(str_detect(fName,paste0("^",prefix))==F) {                  # look for prefix at start of line             
        prefix <- getUserInput(paste(fName, " does not conform to defaults, change prefix from ", prefix," to: "))
        ndPre <- nchar(prefix) + nchar(sep)
        ndCell <- ndPre + nCellDig
      } 
    
      if (nchar(fName) != (ndPre + nCellDig + nTrialDig + nchar(sep) + 4)) {    # either the number of digits in Cell or trial is off
        nCellDig <- getUserInput("how many digits for cell # [default=4]","I")
        if (is.na(nCellDig)) nCellDig <- 4
        nTrialDig <- getUserInput("how many digits for trial # [default=4]","I")
        if (is.na(nTrialDig)) nTrialDig <- 4
        ndPre <- nchar(prefix) + nchar(sep)
        ndCell <- ndPre + nCellDig
        if (nchar(fName) != (ndPre + nCellDig + nTrialDig + nchar(sep) + 4)) {
          print(paste("something's wrong, lets try again, NOTE: file name must end in .ibw"))
          minis()
        }
      } 
      
      cell <- as.integer(substr(fName,(ndPre+1),(ndPre+nCellDig)))
      if (iFil==1) {
        firstCell <- cell
     } else if (cell != firstCell) {
        print(paste("data from multiple cells in folder, only 1 allowed"))
        minis()
      }
      trlStrtChar <- ndCell+nchar(sep)+1
      trlEndChar <- ndCell+nchar(sep)+nTrialDig
      trial[iFil] <- as.integer(substr(fName,trlStrtChar,trlEndChar))
      titl[iFil] <- paste0("C",cell,"T",trial[iFil])
    }
  
    datFils <- data.table(trial=trial, titl=titl, rins= rins, vHold=vHold, nMini=nMini, fil=flist, rms=rms)
    theTrials <- sort(unique(datFils$trial))
    mRes <- list()      # list of datatables of fits and extracted values; one per file, sep list for each cell
    xMins <- list()     # list of arrays of extracted minis; one per file, sep list for each cell
    xMinsf <- list()     # list of arrays of extracted minis; one per file, sep list for each cell
    xMinFit <- list()    # same for fits
    if (bVerbose) print(paste("nFiles=",nFiles, "nTrials=", length(theTrials)))  
    for (t in 1:length(theTrials)) {            # t is counter, trials might not be consecutive
      trl <- theTrials[t]
      print(paste("*******************************************************************************"))
      print(paste("titl=",titl[t]))
      f <- paste0(inFold,"/",datFils$fil[t])
      dList <- getIgorBin(f)
      x <- dList$data
      keyVals <- dList$keyVals
      sealStrt <- as.numeric(keyVals[keys=="SealTestStartX(s)", vals])
      sealLen <- as.numeric(keyVals[keys=="SealTestLength(s)", vals])
      sealAmp <- as.numeric(keyVals[keys=="SealTestAmp(V)", vals])
      if (bMedFilt) {
        xf <- runmed(x, mFiltPts)
      } else {
        xf <- x
      }
#
# Get input resistance and holding current
# NOTE: not median filtering wave first, might be needed for seal test analysis.
# some waves are corrupted with a long string of NaNs, if this is one of those get rid of it.
#
      if (bDoRinVhold) {
#print(paste0("about to call sealRin, sr= ",sr,"keyVals=",keyVals))
        holdAndRin = sealRin(xf, keyVals, sr)
        if (is.nan(holdAndRin[1]) | is.nan(holdAndRin[2]) ) {        # bad data
          print(paste("holdAndRin is bad"))
# for now not moving bad traces, need to test for these values at the end and exclude if averaging
          next
        } else {
          datFils$vHold[t] <- holdAndRin[1]
          datFils$rins[t] <- holdAndRin[2]
        }
      }
#
# Now detect minis
#
      mWaveStrt = 1 + sr*(sealStrt+sealLen+sealPad)
      nTempl <- length(x)-tlen+1  # throw away last bit of data to small to contain a mini template
      dPar$mWaveStrt <- mWaveStrt               # store this in dPar, so we can use it for plotting later
      mf <- xf[(mWaveStrt):(nTempl)]                     # the "mini" containing portion of the file
      m <- x[(mWaveStrt):(nTempl)]                       # the "mini" containing portion of the file
      dlen <- length(mf)
      dPar$dlen <- dlen               # store this in dPar, so we can use it for plotting later
#
# Calculate RMS of m
#
#browser()  
      mSc <- m*aScale
      mnm <- mean(mSc)
      datFils$rms[t] <- mean((mSc-mnm)^2)        # root mean square in pA
#
# scale to pA, change polarity to positive
#
      mf <- mf*aScale*mPol              # NOTE: if need original data it is in X, everthing else is scaled, inverted.
  
      slf <- gsignal::butter(n=4, w=superLow/(sr/2), type="high", plane="z")
    
      bf <- gsignal::butter(n=4, w=cutOff/(sr/2), type="low", plane="z")
    #
# fourth order, freq normalized by half the sampling rate (Nyquist), 'z' means digital rather
# than analog filter; Autoregressive-Moving average seems like default.
### NOTE: vector to be filtered must have mean of 0 (i.e. subtract mean) #######
### ALSO NOTE: even with filtfilt, there was a significant delay between original peak and filtered peak. 
#
      mslf <- filtfilt(slf, mf)         # super low (high pass) filtered version
#
# Detrend: fit a line to entire trace or do low frequency high pass filtering. 
#
      if (detrendMode == "linear") {
        mf <- gsignal::detrend(as.numeric(mf), p=1)
#        print(paste("linear"))
      } else if (detrendMode == "filter") {
        mf <- mslf
#        print(paste("filter"))
      }
      mff <- filtfilt(bf, mf)           # filtfilt performs forward and reverse filtering to minimize delay
# note this filtering done after detrending
    
      tTrace <- (1:dlen)/sr                           # data in tDat is in pA, sec   !!!
      tDat <- data.table(t=tTrace, samp=1:dlen, y=m, yf=mf, yff=mff) # no high pass filter
#fig <- plot_ly(tDat, x=~t, y=~yff, color=I("black"), name='filt', type='scatter', mode='lines')
#fig <- fig %>% add_lines(y=~yf, name='orig ', linetype='dot', alpha=0.2)
#print(fig)
#browser()
      riseLocs <- detectPks(tDat, templ, dPar) # brute force template matching approach
      datFils$nMini[t] <- length(riseLocs)

      if (length(riseLocs) == 0) {
        print(paste("no peaks"))
      } else {
#riseLocs <- riseLocs[1:5]  # diagnostic, just do first 5 peaks        
        sweepRslts <- measurePks(titl[t], tDat, riseLocs, dPar) # return results in list inc DT of vals and array of minis
        mRes[[t]] <- sweepRslts$res1
        xMins[[t]] <- sweepRslts$xMin1
        xMinsf[[t]] <- sweepRslts$xMin1f
        xMinFit[[t]] <- sweepRslts$xFit1
        nMini[t] <- length(riseLocs)       # measurePks should omit if rejecting mini
#
# display--for now display whole trial.
# will change color of points corresponding to rejected and accepted minis
        err <- plotTrl(tDat, titl[t], nMini[t], dPar, mRes[[t]], outdir)
#browser()
      } # end if no minis
    
    } # trial loop
#  
#
# plot rins and vholds
#
    if ((bDoRinVhold) & (nFiles >1)) {
#    datFils <- datFils %>% dplyr::filter(trial <=30)
      p <- ggplot(data = datFils, aes(x=trial, y=rins)) + geom_point() + geom_line(color="black")
      q <- ggplot(data = datFils, aes(x=trial, y=vHold)) + geom_point() + geom_line(color="black")
# NOTE: outdir has terminal backslash
      outf <- paste0(outdir,"Cell",cell,"_rin.pdf")
      outf2 <- paste0(outdir,"Cell",cell,"_vHold.pdf")
      dev.new(pdf)
      pdf(outf, width=10, height=6, pointsize=8)
      plot(p)
    dev.off()
    pdf(outf2, width=10, height=6, pointsize=8)
      plot(q)
      dev.off()
    }

#******************************************************************************************
#
# EXCLUSION STRATEGY: set xclude to code if meets exclusion criteria. Do this on a file by file
#  basis, but keep the minis and values. This way it can be toggled manually in mini display before 
#  calculating final summary values there.
# CHANGING WHEN EXCLUSION OCCURS, so that it is done file by file in measureMini, so that plotting of each file can be colored to reflect it.
# 
    totmins <- sum(nMini)
    if (bVerbose) print(paste("summarizing",totmins,"minis for",length(theTrials),"trials"))
    allMinV <- vector(mode="double", length = mLen*totmins)        #mLen, mini length is defined above
    allMinVf <- vector(mode="double", length = mLen*totmins)        #storing filtered minis also
    allMinFit <- vector(mode="double", length = mLen*totmins)        #and fits
    allMinID <- vector(mode="character", length = mLen*totmins)
    # note, for plotting we need long form where each mini is a group, rather than a column
    pkLoc <- vector(mode="integer", length=totmins)
    TTP <- vector(mode="integer", length=totmins)
    pkAmp <- vector(mode="double", length=totmins)
    rise10 <- vector(mode="integer", length=totmins) 
    rise20 <- vector(mode="integer", length=totmins) 
    rise80 <- vector(mode="integer", length=totmins) 
    rise90 <- vector(mode="integer", length=totmins) 
    mStrt <- vector(mode="integer", length=totmins)
    mTr <- vector(mode="character", length=totmins)   # name minis for trial and pkLoc
    intgrl <- vector(mode="double", length=totmins)
    baseIncpt <- vector(mode="double", length=totmins)
    baseIp <- vector(mode="double", length=totmins)
    baseSlope <- vector(mode="double", length=totmins)
    baseSp <- vector(mode="double", length=totmins)
    baseResid <- vector(mode="double", length=totmins)
    baseMn <- vector(mode="double", length=totmins)
    riseIncpt <- vector(mode="double", length=totmins)
    riseIp <- vector(mode="double", length=totmins)
    riseSlope <- vector(mode="double", length=totmins)
    riseSp <- vector(mode="double", length=totmins)
    riseResid <- vector(mode="double", length=totmins)
    dkAmp <- vector(mode="double", length=totmins)
    dkAmpP <- vector(mode="double", length=totmins)
    dkTau <- vector(mode="double", length=totmins)
    dkTauP <- vector(mode="double", length=totmins)
    dkResid <- vector(mode="double", length=totmins)
    totResid <- vector(mode="double", length=totmins)
    totFit <- vector(mode="double", length=totmins)
    linRatio <- vector(mode="double", length=totmins)
    linIncpt <- vector(mode="double", length=totmins)
    linIp <- vector(mode="double", length=totmins)
    linSlope <- vector(mode="double", length=totmins)
    linSp <- vector(mode="double", length=totmins)
    linResid <- vector(mode="double", length=totmins)
    pAmp <- vector(mode="double", length=totmins)
    dAmp <- vector(mode="double", length=totmins)
    pRise <- vector(mode="double", length=totmins)
    dRise <- vector(mode="double", length=totmins)
    pDK <- vector(mode="double", length=totmins)
    dDK <- vector(mode="double", length=totmins)
    
    xclude <- vector(mode="integer", length=totmins) 
    mCount = 1                            # mini counter for moving values from trials to single vector for cell
#
    for (t in 1:length(theTrials)) {      # t is counter, trials might not be consecutive
      if (nMini[t] > 0) {
        lastMin <- mCount + nMini[t] - 1
        trl[mCount:lastMin] <- theTrials[t]
        pkLoc[mCount:lastMin] <- mRes[[t]]$pkLoc[1:nMini[t]]
        TTP[mCount:lastMin] <- mRes[[t]]$TTP[1:nMini[t]]
        pkAmp[mCount:lastMin] <- mRes[[t]]$pkAmp[1:nMini[t]]
        rise10[mCount:lastMin] <- mRes[[t]]$rise10[1:nMini[t]]
        rise20[mCount:lastMin] <- mRes[[t]]$rise20[1:nMini[t]]
        rise80[mCount:lastMin] <- mRes[[t]]$rise80[1:nMini[t]]
        rise90[mCount:lastMin] <- mRes[[t]]$rise90[1:nMini[t]]
        mStrt[mCount:lastMin] <- mRes[[t]]$mStrts[1:nMini[t]]
        intgrl[mCount:lastMin] <- mRes[[t]]$intgrl[1:nMini[t]]

        baseIncpt[mCount:lastMin] <- mRes[[t]]$baseIncpt[1:nMini[t]]
        baseIp[mCount:lastMin] <- mRes[[t]]$baseIp[1:nMini[t]]
        baseSlope[mCount:lastMin] <- mRes[[t]]$baseSlope[1:nMini[t]]
        baseSp[mCount:lastMin] <- mRes[[t]]$baseSp[1:nMini[t]]
        baseResid[mCount:lastMin] <- mRes[[t]]$baseResid[1:nMini[t]]
        baseMn[mCount:lastMin] <- mRes[[t]]$baseMn[1:nMini[t]]
        riseIncpt[mCount:lastMin] <- mRes[[t]]$riseIncpt[1:nMini[t]]
        riseIp[mCount:lastMin] <- mRes[[t]]$riseIp[1:nMini[t]]
        riseSlope[mCount:lastMin] <- mRes[[t]]$riseSlope[1:nMini[t]]
        riseSp[mCount:lastMin] <- mRes[[t]]$riseSp[1:nMini[t]]
        riseResid[mCount:lastMin] <- mRes[[t]]$riseResid[1:nMini[t]]
        dkAmp[mCount:lastMin] <- mRes[[t]]$dkAmp[1:nMini[t]]
        dkAmpP[mCount:lastMin] <- mRes[[t]]$dkAmpP[1:nMini[t]]
        dkTau[mCount:lastMin] <- mRes[[t]]$dkTau[1:nMini[t]]
        dkTauP[mCount:lastMin] <- mRes[[t]]$dkTauP[1:nMini[t]]
        dkResid[mCount:lastMin] <- mRes[[t]]$dkResid[1:nMini[t]]
        totResid[mCount:lastMin] <- mRes[[t]]$totResid[1:nMini[t]]
        totFit[mCount:lastMin] <- mRes[[t]]$totFit[1:nMini[t]]
        linRatio[mCount:lastMin] <- mRes[[t]]$linRatio[1:nMini[t]]
        linIncpt[mCount:lastMin] <- mRes[[t]]$linIncpt[1:nMini[t]]
        linIp[mCount:lastMin] <- mRes[[t]]$linIp[1:nMini[t]]
        linSlope[mCount:lastMin] <- mRes[[t]]$linSlope[1:nMini[t]]
        linSp[mCount:lastMin] <- mRes[[t]]$linSp[1:nMini[t]]
        linResid[mCount:lastMin] <- mRes[[t]]$linResid[1:nMini[t]]
        pAmp[mCount:lastMin] <- mRes[[t]]$pAmp[1:nMini[t]]
        dAmp[mCount:lastMin] <- mRes[[t]]$dAmp[1:nMini[t]]
        pRise[mCount:lastMin] <- mRes[[t]]$pRise[1:nMini[t]]
        dRise[mCount:lastMin] <- mRes[[t]]$dRise[1:nMini[t]]
        pDK[mCount:lastMin] <- mRes[[t]]$pDK[1:nMini[t]]
        dDK[mCount:lastMin] <- mRes[[t]]$dDK[1:nMini[t]]
        
        xclude[mCount:lastMin] <- mRes[[t]]$xclude[1:nMini[t]]
        
        mTr[mCount:lastMin] <- paste0("T",t,"_",pkLoc[mCount:lastMin])
        mSampl <- 1+ (mCount-1)*mLen                       # sample counter for putting minis consecutively in vector
        allMinV[mSampl:(mLen*lastMin)] <- xMins[[t]]
        allMinVf[mSampl:(mLen*lastMin)] <- xMinsf[[t]]
        allMinFit[mSampl:(mLen*lastMin)] <- xMinFit[[t]]
        mCount <- lastMin + 1
        if (bVerbose) print(paste("mCount=",mCount,"lastMin=",lastMin))
      } # end if nonzero
    } # end for t in theTrials
#
# put into large datatable then separate
# Now separating summary results into two tables (saving to 2 files) with mTr as key. Second table has mostly goodness of fit parameters.
# 
#   
    mTrS <- rep(mTr, each=mLen)            # mini ids repeated for each sample, to make exclusion below easier
    minDT <- data.table(allMinV, allMinVf, allMinFit, mTrS)    # 
    sumDT <- data.table(mTr, trl, pkLoc, TTP, mStrt, rise10, rise90, intgrl, pkAmp, baseSlope, baseMn, riseSlope,dkAmp, dkTau, totResid, linRatio, linResid, pAmp, dAmp, pRise, dRise, pDK, dDK, xclude)
    sumDT2 <- data.table(mTr, rise20, rise80, baseIncpt, baseIp, baseSp, baseResid, riseIncpt, riseIp, riseSp, riseResid, dkAmpP, dkTauP, dkResid, totFit, linIncpt, linIp, linSlope, linSp)
#
# Now separate into included and excluded mini files and results files
#
    Gmins <- sumDT[xclude==0, mTr]         # length of # of good minis, much shorter than mini samples
    Bmins <- sumDT[xclude > 0, mTr]
    GminSamps <- minDT$mTr %in% Gmins      # longer in shorter allowed with vector, but not with data.table
    GminDT <- minDT[GminSamps,]
    BminDT <- minDT[!(GminSamps),]
    GsumDT <- sumDT[mTr %in% Gmins,]
    BsumDT <- sumDT[!(mTr %in% Gmins),]
    GsumDT2 <- sumDT2[mTr %in% Gmins,]
    BsumDT2 <- sumDT2[!(mTr %in% Gmins),]
    #
#
# assign minis to groups, have to assign group number and numInGroup to allMinsDT and to allRes
# switching to a fixed group size of 10--not sure its helpful to look at more than this at once
#
#browser()
    minGrpSize <- 10L
    totminsG <- nrow(GsumDT)
    totminsB <- nrow(BsumDT)
    nGrpsG <- as.integer(ceiling(totminsG / minGrpSize))
    nGrpsB <- as.integer(ceiling(totminsB / minGrpSize))
    grpG  <- vector(mode="integer", length=totminsG)           # don't think these need declaration
    nInGrpG  <- vector(mode="integer", length=totminsG)
    grpB  <- vector(mode="integer", length=totminsB)           # don't think these need declaration
    nInGrpB  <- vector(mode="integer", length=totminsB)
    grpMinG <- vector(mode="character", length = mLen*totminsG)
    grpMinB <- vector(mode="character", length = mLen*totminsB)
    
    lastGrpG <- as.integer((nGrpsG-1)*minGrpSize)             # ending index of last full grp
    if (nGrpsG > 1) grpG[1:lastGrpG] <- rep(1:(nGrpsG-1), each=minGrpSize)    # if there is a full group, repeat each group number except last, minGrpSize times
    grpG[(lastGrpG+1):totminsG] <- nGrpsG
    GsumDT$"grp" <- grpG
    GsumDT2$"grp" <- grpG
    lastGrpB <- as.integer((nGrpsB-1)*minGrpSize)
    if (nGrpsB > 1) grpB[1:lastGrpB] <- rep(1:(nGrpsB-1), each=minGrpSize)    # repeat each group number except last, minGrpSize times
    grpB[(lastGrpB+1):totminsB] <- nGrpsB
    BsumDT$"grp" <- grpB
    BsumDT2$"grp" <- grpB
    
# Add columns to summary DTs: nInGrp is which number the mini is within group, grp is group #, 
    
    nInGrpG <- rep(1:(minGrpSize), length.out=totminsG)       # this is for summary files
    nInGrpB <- rep(1:(minGrpSize), length.out=totminsB)       # 
    GsumDT$"nInGrp" <- nInGrpG
    BsumDT$"nInGrp" <- nInGrpB
    GsumDT2$"nInGrp" <- nInGrpG
    BsumDT2$"nInGrp" <- nInGrpB
    
# Add columns to mini DTs: numInGrp is which number the mini is within group, grpMin is group #,
# samps is sample number, allMinID is the mTr which encodes the sample # of the start and the trial
# number. For all these columns, need to repeat values for each sample in each mini
    
    grpMinG <- rep(grpG, each = mLen)
    grpMinB <- rep(grpB, each = mLen)
    numsIn1Grp <- rep(1:minGrpSize, each=mLen)             
    numInGrpG <- rep(numsIn1Grp, length.out=nrow(GminDT)) # 
    numInGrpB <- rep(numsIn1Grp, length.out=nrow(BminDT)) #

    GminDT$"grpMin" <- grpMinG 
    BminDT$"grpMin" <- grpMinB 
    GminDT$"numInGrp" <- numInGrpG
    BminDT$"numInGrp" <- numInGrpB
    
    sampsG <- rep.int(c(1:mLen),totminsG)
    sampsB <- rep.int(c(1:mLen),totminsB)
    GminDT$"samps" <- sampsG
    BminDT$"samps" <- sampsB
#
#  may want to include cell name in summary file and mini file names. 
#
#browser()
    cellAmp[iCell] <- round(mean(GsumDT$pkAmp),1)
    cellGdMins[iCell] <- totminsG
print(paste("***FOR CELL",cellStr[iCell], "#good minis=",cellGdMins[iCell], "avg amp=",cellAmp[iCell]))

    outf <- paste0(outdir,cellStr[iCell],"_datFils.txt")
    fwrite(datFils, outf ,sep = "\t")
    outf <- paste0(outdir,cellStr[iCell],"_dPar.txt")              # list of detection parameters, can be written, read like DT
    fwrite(dPar, outf ,sep = "\t")


    summOut <- paste0(outdir,cellStr[iCell],"_filtSumG.txt")
    fwrite(GsumDT, summOut ,sep = "\t")
    summOut <- paste0(outdir,cellStr[iCell],"_filtSumB.txt")
    fwrite(BsumDT, summOut ,sep = "\t")
    summOut <- paste0(outdir,cellStr[iCell],"_filtSumG2.txt")
    fwrite(GsumDT2, summOut ,sep = "\t")
    summOut <- paste0(outdir,cellStr[iCell],"_filtSumB2.txt")
    fwrite(BsumDT2, summOut ,sep = "\t")

    minOut <- paste0(outdir,cellStr[iCell],"_filtMinisG.txt")
    fwrite(GminDT, minOut, sep = "\t")
    minOut <- paste0(outdir,cellStr[iCell],"_filtMinisB.txt")
    fwrite(BminDT, minOut, sep = "\t")
  
  } # next cell
  cellSumm <- data.table(cellStr, cellAmp, cellGdMins)
  fwrite(cellSumm, paste0(conditionFold,"/cell_summary.txt"), sep = "\t")

} # end minis function
#**********************************************************************
#
# miniStat (vMin, strt, pkLoc, tauDK, nTau)
# calculates probabilities that 1) rise + decay phase are > baseline phases, 2) that mean derivatives of rise are > mean derivatives of baseline, and that
# mean derivatives of decay are < mean derivatives of baseline. For each contrast, perform Welch's t-test and compute Cohen's D.
#
# vMin is a real vector containing the full mini; strt and pkLoc are indices into that vector of the first point after the initial baseline and the point of the peak.
# the end of the mini is defined by the tauDK and nTau (number of taus) to include past the peakLoc. This should leave additional time for a second baseline period.
# For example, a monoexponential function decays to 0.37^3 = 0.05 in 3 time constants and we are extracting ~ 5 taus (of the template anyway)
#
# Returns list containing dataFrame of results and error code. DataFrame has 3 rows for 3 comparisons, columns for comparison type (Amp, Rise, Decay), for p and for D. 
#
#**********************************************************************
  miniStat <- function (vMin, strt, pkLoc, tauDK, nTau) {
  
    minNd <- pkLoc + round((nTau*tauDK), 0)                    # end of the non-baseline portion of the mini
    mLen <- length(vMin)
    if ((mLen < minNd) | (mLen < (pkLoc + round((nTau*2.3),0)))) {      # if < time requested or < time needed to decay to 10% of amplitude
      results <- list("pd"=NULL, "err"=-1)
      return (results) 
    }
    if (((strt < 1) | (minNd < 1) | ((pkLoc-strt)<4))) {  # need positive non zero vals for strt,minNd and at least 3 rise pts= 2 rise slope pts
      results <- list("pd"=NULL, "err"=-1)
      return (results) 
    }
    
    vBase <- c(vMin[1:(strt-1)],vMin[minNd:mLen])               # add together the pre and post mini baselines
    vRise <- vMin[strt:pkLoc]
    vDK <- vMin[(pkLoc+1):minNd]
    
    dRise <- shift(vRise, type="lead")[1:(length(vRise)-1)] - vRise[1:(length(vRise)-1)]      # implementation of sample by sample derivative--defined for all but first sample
    dBase <- shift(vBase, type="lead")[1:(length(vBase)-1)] - vBase[1:(length(vBase)-1)]
    dDK <- shift(vDK, type="lead")[1:(length(vDK)-1)] - vDK[1:(length(vDK)-1)]
    
# test means of rise + decay are > baseline
    tAmp <- try(t.test(c(vRise,vDK), vBase, alternative="g"))
    if(inherits(tAmp, "try-error") ) {
      print(paste("error amp t-test"))      
      results <- list("pd"=NULL, "err"=-1)
      return (results) 
    } else {
      sdAmp <- sqrt((var(c(vRise,vDK)) + var(vBase))/2)
      tidAmp <- tidy(tAmp)
      fxAmp <- abs(tidAmp$estimate2[1] - tidAmp$estimate1[1])/sdAmp
    }
# test mean rise derivs are > baseline derivs
    tdRise <- try(t.test(dRise, dBase, alternative="g"))
    if(inherits(tdRise, "try-error") ) {
      print(paste("error rise t-test"))      
      results <- list("pd"=NULL, "err"=-1)
      return (results) 
    } else {
      sdRise <- sqrt((var(dRise) + var(dBase))/2)
      tidRise <- tidy(tdRise)
      fxRise <- abs(tidRise$estimate2[1] - tidRise$estimate1[1])/sdRise
    }
# test mean decay derivs are < baseline derivs
    
    tdDK <- try(t.test(dDK, dBase, alternative="l"))
    if(inherits(tdDK, "try-error") ) {
      print(paste("error DK t-test"))      
      results <- list("pd"=NULL, "err"=-1)
      return (results) 
    } else {
      sdDK <- sqrt((var(dDK) + var(dBase))/2)
      tidDK <- tidy(tdDK)
      fxDK <- abs(tidDK$estimate2[1] - tidDK$estimate1[1])/sdDK
    }
    comparType <- c("Amp","Rise_Slope","Decay_Slope")
    pWelchT <- signif(c(tidAmp$p.value[1], tidRise$p.value[1], tidDK$p.value[1]), 2)
    cohensD <- round(c(fxAmp, fxRise, fxDK),2)
    pd <- tibble(comparType, pWelchT, cohensD)
    
    results <- list("pd" = pd, "err"=0)
#browser()
    return (results) 
  }
#**********************************************************************
minis() # the main loop is now a function so I can recall it
