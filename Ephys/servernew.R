


shinyServer(function(input, output, session) {
  
#  allMinsDT$grp <- as.integer(allMinsDT$grp)
  minGrpSize <- 10L 
  cols <- brewer.pal(11,"Spectral")
  pal <- colorRampPalette(cols)
  hexcols <- pal(minGrpSize)
#  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  
  observeEvent (ignoreInit=T, c(input$group, input$checkGroup, input$bFit, input$rb, input$gOrb),
    { minGroup <- input$group
#browser()
      if (is.null(input$checkGroup)) {
        minNums <- "1"
      } else {
        minNums <- input$checkGroup
      }
      if (input$gOrb == "Good") {
        theGroup <- gMinDT[grpMin == minGroup,]   # first select all the data for this group of minis
        theSum <- gSumDT
      } else {
        theGroup <- bMinDT[grpMin == minGroup,]   # first select all the data for this group of minis
        theSum <- bSumDT
      }
      v <- theGroup$numInGrp %in% minNums    # select vector that match desired mini nums
      theData <- theGroup[v,]
      titl <- paste(theData[1, mTrS], collapse=" ")
      nPlots <- length(minNums)
      mLen <- max(theData$samps)
      for (i in 1:nPlots) {
        strt <- 1 + (i-1)*mLen
        nd <- strt+mLen-1
        theMini <- theData[(strt):(nd),]
        if (i==1) {
          if (input$rb=="filtered") {
            p <- plot_ly(theMini, x = ~samps) %>%
              add_trace(type="scatter", y = ~allMinVf, mode = 'lines')
          } else {
            p <- plot_ly(theMini, x = ~samps) %>%
              add_trace(type="scatter", y = ~allMinV, mode = 'lines')
          }
        } else {
          if (input$rb=="filtered") {
            p <- p %>% add_trace(data = theMini, type="scatter", x = ~samps, y = ~allMinVf, mode = 'lines')
          } else {
            p <- p %>% add_trace(data = theMini, type="scatter", x = ~samps, y = ~allMinV, mode = 'lines')
          }
        }
        if (input$bFit) {
          p <- p %>% add_trace(data=theMini,x = ~samps, y = ~allMinFit,
            type='scatter', mode='markers')
        }
      }
      p <- layout(p,showlegend=F, title = titl)   # note, not bothering with titles for multiple minis
      output$mPlot <- renderPlotly({p})
      output$mTable <- renderDT(datatable(req(theSum), 
        extensions = 'Buttons', options = list(dom='Bfrtip',
        buttons=c('copy', 'csv', 'excel', 'print', 'pdf'), 
        pageLength = minGrpSize, displayStart=(minGroup-1)*minGrpSize,
        columnDefs = list(list(visible=T, targets=c(1:25)))),
        editable = list(target = "cell", disable = list(columns = c(1:10)))))
        
        
    })

  
  
  observeEvent(input$reload, {
#    shinyDirChoose(input, "directory", roots=c(wd='.'), filetypes=c(''),
#                   defaultPath='', defaultRoot='wd', allowDirCreate = FALSE)
#browser()
#    resFold <- getLastRes(input$directory)
#    pth <- paste0(resFold,"/",basename(inFold),"_")
#    gSumDT <- fread(paste0(pth,"filtSumG.txt"))
#    bSumDT <- fread(paste0(pth,"filtSumB.txt"))
#    gMinDT <- fread(paste0(pth,"filtMinisG.txt"))
#    bMinDT <- fread(paste0(pth,"filtMinisB.txt"))
  
    })
  


    
    
  
  
  
  
  
  })


