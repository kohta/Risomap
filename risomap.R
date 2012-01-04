################################################################################
##　risomap interface for R
##  created by Kohta Ishikawa (quantumcorgi _at_ gmail.com)
##  date       Jan 3 2012
##  2012 Kohta Ishikawa all rights reserved.
##
##  There might be some bugs in this program.
##  You could tell me if you find it.
################################################################################

dyn.load("risomap.so")

#x: data matrix
#k: edge number per data point used in creating knn graph
risomap <- function(x,k){
  if(class(x) != "matrix"){
    x <- as.matrix(x)
  }
  .Call("risomap",x,k)
}

