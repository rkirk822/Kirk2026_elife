# Current version Mar 31
# multiple bug fixes, added cell and condition labels to summary files
# removed capacity to remake filtered good/bad minis files.
#

# Current version: Feb 11, 2023
# fixed inclusion of poor decay fits
# added default for nSamp


library(data.table)
library(tidyverse)
library(plotly)
library(patchwork)


#************************************************************************
# Produce an average mini from a saved results file. This version produces an average mini from all
# of the non-excluded minis of all the cells for each condition (multi cell) 
# Will also produce cumulative histograms of all minis (single cell) 
# for amplitude, rise time, decay time and inter-mini-interval

# nSamp = 0 disables sampling or uses nSamp = minimum common number
#
# Cell folders without results folders throws an error.  
#
#************************************************************************
avgMini <- function(nSamp=0) {
  
# April 23: From now on, assume dPar has been written to results folder, recording detection/extraction parameters.
# Hence will not plan on additional filtering here for maxTTP, maxDkP. If more filtering needed, should probably
# be implemented in shiny routines where individual minis can be visualized and parameters are listed.
# Should no longer need to align minis, since this is done at extraction--removing this routine.
  
  prefix <- "Cell"; nCells <- 0
  cellOrCondition <- getUserInput(paste("All cells for one condition [return] or single cell mode [any other key]"))        # default varType is string, so not needed
  if (cellOrCondition=="") {
    bMult <- T
    inFold <- rstudioapi::selectDirectory(caption = "Select Condition folder", label = "Select")
    if (is.null(inFold)) break
    cellFolds <- list.dirs(path=inFold, recursive = F)
    nCells = sum(str_detect(cellFolds,regex(prefix,ignore_case=T)))
    if (nCells==0)  {
      print(paste("No enclosed cell folders found"))
      break
    } else {
      cellFolds <- cellFolds[which(str_detect(cellFolds,regex(prefix,ignore_case=T)))]     # gets rid of any additional folders, not cell folders
      cellNums <- str_split_i(cellFolds, regex(prefix,ignore_case=T),i=-1)
    }
    condStr <- basename(inFold)
  } else {            # single cell mode
    bMult <- F
    nCells = 1
    cellFolds <- rstudioapi::selectDirectory(caption = "Select Cell folder", label = "Select")
    if ((is.null(cellFolds)) | str_detect(cellFolds,regex(prefix,ignore_case=T)) == F) break
    condStr <- basename(cellFolds)
  }
  useAll <- F
  if (nSamp==0) {
    useSamp <- getUserInput("Use max common number of minis?[return] or use all minis[any other key]")
    if (useSamp!="") useAll <- T
  }
#
# set up vectors for average amp and interval and # minis to sample
#
  cellMnAmp <- vector(mode="double", length=nCells)
  cellMnISI <- vector(mode="double", length=nCells)
  cellMnRMS <- vector(mode="double", length=nCells)
  nSamps <- vector(mode="integer", length=nCells)
#
# if this is multicell, find number of minis in each trial.
#
  if (bMult==T) {
    nGminis <- vector(mode="integer", length=nCells)
    for (iCell in 1:nCells) {
      lastRes <- getLastRes(cellFolds[iCell])
      cellStr <- basename(cellFolds[iCell])
      pth <- paste0(lastRes,"/", cellStr,"_filtSumG.txt")
      GsumDT <- fread(file=pth)
      nGminis[iCell] <- nrow(GsumDT)
print(paste("#good minis=",nGminis[iCell], "for file ", cellStr,"_filtSumG.txt"))
      dPar <- fread(paste0(lastRes,"/", cellStr,"_dPar.txt"))        # get detection/extraction parameters
    }
    if (useAll == F) {
      nSamps[1:(nCells)] <- min(nGminis)
      print(paste("using nSamps=",nSamps[1]," min from cell",cellFolds[which.min(nGminis[iCell])]))
      nSampRows <- nCells*nSamps[1]
    } else {                          # using all minis, store sum in nSampRows
      nSamps <- nGminis                # and use each number of minis as the "sample"
      nSampRows <- sum(nGminis)
    }
  } else {                            # just one cell, but need to get the number of minis and cell number string
    lastRes <- getLastRes(cellFolds[1])
    cellNums <- str_split_i(cellFolds, regex(prefix,ignore_case=T),i=-1)
    GsumDT <- fread(paste0(lastRes,"/",basename(cellFolds[1]),"_filtSumG.txt"))
    nSampRows <- nrow(GsumDT)
    nSamps <- nrow(GsumDT)
  }
#
#  Now, loop through one or more cells
#
#browser()
  for (iCell in 1:nCells) {
    lastRes <- getLastRes(cellFolds[iCell])
    cellFname <- basename(cellFolds[iCell])
    pth <- paste0(lastRes,"/",cellFname,"_filtMinisG.txt")
    GminDT <- fread(file=pth)
    pth <- paste0(lastRes,"/",cellFname,"_filtMinisB.txt")
    BminDT <- fread(file=pth)
    pth <- paste0(lastRes,"/",cellFname,"_filtSumG.txt")
    GsumDT <- fread(file=pth)
    pth <- paste0(lastRes,"/",cellFname,"_filtSumB.txt")
    BsumDT <- fread(file=pth)
    pth <- paste0(lastRes,"/",cellFname,"_datFils.txt")
    datFils <- fread(file=pth)
    cellMnRMS[iCell] <- round(mean(datFils$rms),1)
    if (iCell==1) {
      aggDatFils <- datFils
    } else {
      aggDatFils <- rbind(aggDatFils, datFils)
    }
#
# Need to turn pkLoc (in samples) into ISI in msec. pkLoc is in samples since start of file, but files are
# put together in summary, so pkLoc will jump back to beginning of 10s period with every file. This results
# in negative ISIs. Replace these with the pkLoc value (ISI since beginning of trial)
#
    GsumDT$"ISI" <- (GsumDT$pkLoc - shift(GsumDT$pkLoc, fill=0, type="lag"))/10
    GsumDT[ISI<0, ISI := (pkLoc/10)]
    GsumDT$"cellNum" <- cellNums[iCell]         # add column to keep track of which mini is from which cell
#browser()
#
# calculate cell averages
#
    cellMnAmp[iCell] <- round(mean(GsumDT$pkAmp),1)
    cellMnISI[iCell] <- round(mean(GsumDT$ISI),1)
    mLen <- max(GminDT$samps)    # the length of minis in the input file
    print(paste("# good minis=",nrow(GsumDT)))
    print(paste("input mini length=",mLen,"nSamps=",nSamps[iCell]))
    minSamps <- vector(mode = "integer", length=nSamps[iCell]*mLen)
    rows2samp <- 1:(nrow(GsumDT))
    samps <- sample(rows2samp, nSamps[iCell])  # note: excluding first trials from sampling
    for (iSamp in 1:nSamps[iCell]) {              # compute row numbers in long format mini file
# NOTE, we are using unpadded mini length, not yet padded!!!
      minSamps[(1+(iSamp-1)*mLen):(iSamp*mLen)] <- (1:mLen)+(samps[iSamp]-1)*mLen
    }
#
# If this is the first cell (and we are sampling) preallocate DTs to hold sampled parameters and minis
# To preallocate we first get the column names and types for GminDT and GsumDT. 
#
    if (bMult == F) {              # single cell, no sampling
      aggSum <- GsumDT
      aggMinIn <- GminDT
      print(paste("folder: ",cellFolds[iCell], "cell # ", iCell))
      titl <- inFold
    } else {                       # multi-cell,
#
# we are sampling, either with a fixed number, or using all of the minis
#
      rows2samp <- 1:(nrow(GsumDT))
      samps <- sample(rows2samp, nSamps[iCell])  # 
      titl <- cellFolds[iCell]
      # If this is the first cell, set up the aggregation data.tables
      if (iCell==1) {
        sumCols <- sapply(GsumDT, class)    # returns a named character vector of column types, names are column names
        minCols <- sapply(GminDT, class)
        aggSum <- data.table()
        aggMinIn <- data.table()
        for (icol in 1:length(sumCols)) {   # recreate the vectors, now with the length we will need                    
          v <- vector(mode = sumCols[icol], length=nSampRows)  # make a vector of the type of the corresponding column in GsumDT
          aggSum[,(names(sumCols)[icol])] <-v   # assign v to aggSum column, with the name of the corresponding column in GsumDT
        }
        for (icol in 1:length(minCols)) {   # recreate the vectors, now with the length we will need                    
          v <- vector(mode = minCols[icol], length=(nSampRows*mLen))
          aggMinIn[,(names(minCols)[icol])] <-v   # assign v to aggMin column, with the name of the corresponding column in GminDT
        }
        #browser()
        aggSum[1:(nSamps[iCell]),] <- GsumDT[samps,]
        aggMinIn[1:(mLen*nSamps[iCell]),] <- GminDT[minSamps,]
      } else {                    # this is not cell 1, and we ARE doing multiple cells
#        browser()
#
# each sample of the summary DT takes nSamp rows; each sample of minis takes nSamp minis times mLen values per mini
#
        if (useAll) {                              # nSamps is number of minis, differs from cell to cell
          strtS <- 1 + sum(nSamps[1:(iCell-1)])    # start the next summary at the row following the sum of the minis in prior cells
          strtM <- 1 + sum(nSamps[1:(iCell-1)])*mLen

        } else {                                   # not using all, sampling same number from each
          strtS <- 1+(iCell-1)*nSamps[iCell]
          strtM <- 1+((iCell-1)*mLen*nSamps[[iCell]])
        }
        ndS <- strtS + nSamps[iCell] - 1              # end 1 row less than # minis in this Cell later
        ndM <- strtM + nSamps[iCell]*mLen - 1         # end 1 less than # minis X length of each mini later
        aggSum[(strtS):(ndS),] <- GsumDT[samps,]
        aggMinIn[(strtM):(ndM),] <- GminDT[minSamps,]
      }   # end if icell==1
    }     # end if bMult
  }       # end iCell for loop
  #
  # Write aggSum file (for now do not bother with aggregated bad minis). Will need to examine this for goodness of fit etc to 
  # improve efficiency of rejecting minis
  #
  outf <- paste0(inFold,"/",condStr,"_aggGminSum.txt")        # aggregated summaries of all good minis in condition
  fwrite(aggSum, outf, sep="\t")
  outf <- paste0(inFold,"/",condStr,"_aggDatFils.txt")        # aggregated summaries of data files
  fwrite(aggDatFils, outf, sep="\t")
  
  
  #
  # Create grand average mini and cum histos
  # NOTE: not creating avg minis for each cell, but grand avg is avg of one cell in single cell mode
  #
#  browser()
#  midRise <- round(((aggSum$rise20 + aggSum$rise80)/2),0)
#  offs <- aggSum$TTP - midRise
#  offs[offs<0] <- 0
#  minAlign <- alignMinis(aggMinIn, offs, mLenO, 0)
  #browser()
  minAvg <- aggMinIn %>% group_by(samps) %>% summarise(mn = mean(allMinV), mnf = mean(allMinVf))
  mFig <- plot_ly(minAvg, x = ~samps, y = ~mn, type='scatter', mode='lines')
  l <- layout(mFig, title = titl)
  print(mFig)
  outf <- paste0(inFold,"/",condStr,"_avgMini.html")
  # htmlwidgets fails if not connected to the internet
  htmlwidgets::saveWidget(as_widget(partial_bundle(l)), outf, selfcontained = T)
  #
  # saveWidget is supposed to save a self contained file without also saving a folder of other stuff.
  # the index file is self contained, but the extra folder is still created, so just get rid of it.
  #
  outfold <- str_replace(outf,".html","_files")
  unlink(outfold, recursive = T)
  outfold <- str_replace(outfold,"_files",".txt")
  fwrite(minAvg,outfold)

  cellMeans <- data.table(cellMnISI,cellMnAmp,cellMnRMS, cellFolds)
  outf <- paste0(inFold,"/",condStr,"_cellMeans.txt")
  fwrite(cellMeans,outf)
  print(paste("mean Amp=",mean(cellMnAmp)))
  print(paste("mean ISI=",mean(cellMnISI)))
#browser()  
  sortAmps <- sort(aggSum$pkAmp)
  cumAmps <- cumsum(sortAmps)
  sortISI <- sort(round(aggSum$pkLoc/10))
  cumISI <- cumsum(sortISI)    # rounding to nearest msec
#  cumHalfRise <- cumsum(order(round(0.05*(aggSum$rise80+aggSum$rise20),1)))   #  *.05 = averaging and dividing by 10 to msec
  sortIntgrl <- sort(aggSum$intgrl)
  cumIntgrl <- cumsum(sortIntgrl)
  aggSum[dkTau > dPar$maxTDk, dkTau:=0.001]                         # filtering out long tau, code is 0.001
  sortDkTau <- sort(round(aggSum$dkTau,1))
  cumDkTau <- cumsum(sortDkTau)
  iMini <- 1:(nSampRows)
  cums <- data.table(iMini,sortAmps, cumAmps, sortISI, cumISI, sortIntgrl, cumIntgrl, sortDkTau, cumDkTau)
  outf <- paste0(inFold,"/",condStr,"_cumHistos.txt")
  fwrite(cums,outf)
  p1 <- ggplot(cums, aes(x=iMini, y=sortAmps)) + geom_line()
  p2 <- ggplot(cums, aes(x=iMini, y=sortISI)) + geom_line()
  p3 <- ggplot(cums, aes(x=iMini, y=sortIntgrl)) + geom_line()
  p4 <- ggplot(cums, aes(x=iMini, y=sortDkTau)) + geom_line()
  
  outf <- paste0(inFold,"/",condStr,"_cumAmpISI.pdf")
  dev.new(pdf)
  pdf(outf, width=10, height=6, pointsize=8)
  patch <- wrap_plots(p1,p2) & theme_minimal()
  print(patch)
  dev.off()
  outf <- paste0(inFold,"/",condStr, "_cumStuff.pdf")
  dev.new(pdf)
  pdf(outf, width=10, height=6, pointsize=8)
  patch2 <- p3 + p4 & theme_minimal()
  print(patch2)
  dev.off()
}

