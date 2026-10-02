library(shiny)
library(data.table)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(DT)
library(shinyFiles)

inFold <- rstudioapi::selectDirectory(caption = "Select Cell Folder", label = "Select Cell Folder")
resSel <- getUserInput(paste("Most recent results folder [return] or select earlier results folder [any other key]"))
if (resSel=="") {
  resFold <- getLastRes(inFold)
} else {
  resFold <- rstudioapi::selectDirectory(caption = "Select Results Folder", label = "Select Results Folder", path=inFold)
  if (is.null(resFold)) break
}
pth <- paste0(resFold,"/",basename(inFold),"_")
gSumDT <- fread(paste0(pth,"filtSumG.txt"))
if (nrow(gSumDT)==0) {
  print(paste("NO GOOD MINIS"))
  stopApp()
}
bSumDT <- fread(paste0(pth,"filtSumB.txt"))
gMinDT <- fread(paste0(pth,"filtMinisG.txt"))
bMinDT <- fread(paste0(pth,"filtMinisB.txt"))

#
# NOTE: this mini viewer is designed to work with data analyzed by mini4.R and its derivatives
#