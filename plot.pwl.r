# Function that plots the piecewise equations along with the dataset

plot.pwl <- function(pwl, data){
  if(class(pwl)!="pwl") stop("Need to provide the piecewise linear equations")
  
  x <- data[,1]
  y <- data[,2]
  
  #plot the data
  plot(x,y, ylim=c(min(y), max(y)), xlim=c(min(x), max(x)))
  
  BP <- pwl$BreakPoints
  noOfBP <- length(pwl$BreakPoints)
  
  i <- 1
  repeat{
    if(i == 1){
      # plotting the first segment
      curve(pwl$coeffs[1,i] + pwl$coeffs[2,i]*x, add=T, from=0, to=BP[i], col="red")
      abline(v=BP[i], lty=3)
    }else if(i == noOfBP+1){
      # plotting the last segment
      curve(pwl$coeffs[1,i] + pwl$coeffs[2,i]*x, add=T, from=BP[noOfBP], to=max(x), col="red")
    }else{
      # plotting all the other segments
      curve(pwl$coeffs[1,i] + pwl$coeffs[2,i]*x, add=T, from=BP[i-1], to=BP[i], col="red")
      abline(v=BP[i], lty=3)
    }
    i <- i+1
    if(i>noOfBP+1) break()
  }
}