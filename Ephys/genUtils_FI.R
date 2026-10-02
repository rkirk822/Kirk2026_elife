#
# genUtils.R  some general purpose utilities. 
# to use these in your code use: source("genUtils.R") if this file is in the same directory as your code, or
# source("directory_where_you_put_it/genUtils.R") if you put it in "/directory_where_you_put_it"
#

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
