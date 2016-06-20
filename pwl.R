# Find piecewise linear approximation function for the given data and number of breakpoints (bp)
# Length (l) specifies the minimum lenght of the segments (i.e. length between breakpoints)
# class "pwl" consists of the locations of the breakpoints, coefficients of the equations,
# fitted values, residuals, and mse.
# The function returns a list of pwl as result.


pwl <- function(data, noOfBP, l, error, maxBP, ...){

  # check if l is given
  if(missing(l)) stop("Please specify mininum length between break points:\n length is to be at least 1", call. = FALSE)
  
  # check if either no of BP or error is given
  if(missing(noOfBP)&& missing(error) && missing(maxBP)) stop("Please specify either number of desired breakpoints or error", call. = FALSE)
  
  # first sort the given dataset
  data <- data[sort.list(data[,1]),]
  
  # check if the length and no of BP given will fit into the given data
  size <- nrow(data)
  maxBP.allowed <- (size/l) -1
  if(noOfBP > maxBP.allowed) stop("The data set is not big enought to fit the number of breakpoints given.\n Either lower the number of desired break points or distance between break points", call. = FALSE)
  
  result <- list(minssr=0, minmse=0, BP=c())

  #Use this if noOfBP is given..
  
  # If there is only one breakpoint, we don't need to calculate the MSE matrix
  # MSE matrix is used when there are more than 1 BP and/or when error is given
  # and no. of BP is unknown
  
  if(noOfBP == 1){
    result <- findoneBP(data, l)
    result$minmse <- result$minssr/nrow(data)

  }else{
    
    ssrMatrix <- calculateSSRMatrix(data)
    print("SSRMatrix done!")
    result <- findBP(ssrMatrix, noOfBP, l, 1, nrow(data))
    print(paste0("BP found, BP = ", result$BP))
    
    BP <- result$BP
    result$minmse <- result$minssr/nrow(data)
    result$BP <- data[BP,1]
  }

  piecewise <- getequations(data, result$BP)
  piecewise$mse <- result$minmse
  
  class(piecewise) <- "pwl"
  
  #allpwl <- list(piecewise)
  piecewise
  #allpwl
}



