#
# genUtils.R  some general purpose utilities. 
# to use these in your code use: source("genUtils.R") if this file is in the same directory as your code, or
# source("directory_where_you_put_it/genUtils.R") if you put it in "/directory_where_you_put_it"
#
# 2-11-23 added getInput()

#************************************************************************


expDecay1 <- function(t, amp, tau) {
  exp <- amp*(exp(-t/tau))   
}
#************************************************************************

expDecay2 <- function(t, amp, fr, t1, t2) {
  #
  # from https://martinlab.chem.umass.edu/r-fitting-data/
  #
  exp2 <- amp*((fr*exp(-t/t1)) + ((1-fr)*exp(-t/t2)))  
}
#************************************************************************
#fitLine <- function (xy,)

#************************************************************************
remOutlieIQR <- function (v) {
#find Q1, Q3, and interquartile range for values in column A
Q1 <- quantile(v, .25)
Q3 <- quantile(v, .75)
IQR <- Q3 - Q1

#only keep values that have values within 1.5*IQR of Q1 and Q3
r <- v[v> (Q1 - 1.5*IQR) & v< (Q3 + 1.5*IQR)]
return(r)
}
#************************************************************************
# returns a string containing the current day and time, but without specifying the timezone
nowNoTimeZone <- function() {
  return(str_sub(now(), 1, str_length(now())-3))
}
#************************************************************************
# Get user input and do conversion if needed
# varTypes are "S" string--no conversion, "I" integer, "N" number
getUserInput <- function(prompt, varType = "S") {
  theVar <- readline(prompt=prompt)
  if (varType=="I") {
    theVar <- as.integer(theVar)
  } else if (varType=="N") {
    theVar <- as.numeric(theVar)
  }
  return(theVar)
}
#************************************************************************
