createRoll <- function(n){
  theta <- seq(0,2*pi,2*pi/(n-1))
  r <- 5*(1-0.5*theta/(2*pi))
  x <- r*cos(theta)
  y <- r*sin(theta)
  d <- matrix(c(x,y),ncol=2)
}

createRoll3d <- function(){
  createOneLine <- function(y,n){
    tmp <- createRoll(n)
    tmp <- cbind(tmp,rep(y,n),1:n)
    tmpz <- tmp[,2]
    tmp[,2] <- tmp[,3]
    tmp[,3] <- tmpz
    return(tmp)
  }
  rolldata <- lapply(1:10,function(yn){createOneLine(yn/10,30)})
  result <- rolldata[[1]]
  for(i in 1:10){
    result <- rbind(result,rolldata[[i]])
  }
  result[order(result[,4]),-4]
}
