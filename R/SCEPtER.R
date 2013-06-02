errorObs <- function(sigma, STAR) {

  if(length(sigma) != 5)
    stop("uncorrect number of elements in the error vector")
  
  mysigma <- sigma
  # sigma[4,5] are expressed as % ...
  mysigma[4] <- mysigma[4]*STAR[4]
  mysigma[5] <- mysigma[5]*STAR[5]

  # sample the observation from multivariate normal distribution
  # assuming diagonal covariance matrix
  STARp <- MASS::mvrnorm(1, STAR[1:5], diag(mysigma^2))

  # attach the information about the object...
  STARp <- c(STARp, STAR[c(6:9, 3)])

  return(STARp)
}


estimatePert <- function(data, STAR, sigma, thr, sel, Nrun=10000) {

  if(!is.matrix(data))
    stop("the recovery grid must be of \"matrix\" class")
  if(!is.vector(STAR))
    stop("the star target must be a vector")
  if(!is.vector(sigma))
    stop("the sigma object must be a vector")
  if(dim(data)[2] != 9)
    stop("uncorrect number of columns in the data grid")
  if(length(sigma) != 5)
    stop("uncorrect number of elements in the uncertainties vector")
  if(length(sel) != 5)
    stop("uncorrect number of elements in the selection vector")
  
  mysigma <- sigma
  mysigma[4] <- mysigma[4]*STAR[4]
  mysigma[5] <- mysigma[5]*STAR[5]
  
  STARp <- MASS::mvrnorm(Nrun, STAR[1:5], diag(mysigma^2))
  
  v <- NULL
  
  for(i in 1:Nrun) {

    res <- .Call("lik2", data, STARp[i,], sigma, thr, sel)
    
    if( !is.null(res) ) {
      v <- rbind(v, res)
    }
  }
  v <- as.data.frame(v, row.names=1:Nrun)
  colnames(v) <- c("M", "R", "num")

  return(v)
}

estimateError <- function(data, STAR, sigma, thr, sel) {

  # main function: estimate M and R assuming error on observations

  # sanity checks ...
  if(!is.matrix(data))
    stop("the recovery grid must be of \"matrix\" class")
  if(dim(data)[2] != 9)
    stop("uncorrect number of columns in the data grid")
  if(length(sigma) != 5)
    stop("uncorrect number of elements in the error vector")
  if(length(sel) != 5)
    stop("uncorrect number of elements in the selector vector")

  if(!is.matrix(STAR)) 
    STAR <- matrix(STAR, nrow = 1)
    
# likelihood computation and mean of the best points
  dim <- nrow(STAR)

  res <- NULL
  for(i in 1:dim) {

    # generate a "perturbed" object
    STARp <- errorObs(sigma, STAR[i,])

    # ML estimation of Mass and Radius
    val <- .Call("lik2", data, STARp, sigma, thr, sel)

    if(!is.null(val)) {
      val <- c(val, STAR[i, c(6:9, 3)])
      res <- rbind(res, val)
    }
  }
  res <- as.data.frame(res, row.names=1:nrow(res))
  colnames(res) <- c("M", "R", "num", "trueM", "trueR", "logAge", "pcAge", "FeH")
  res$errorM <- (res$M - res$trueM)/res$trueM
  res$errorR <- (res$R - res$trueR)/res$trueR

  return(res)
} 



sampleStar <- function(n, grid) {

  if(!is.matrix(grid))
    stop("the recovery grid must be of \"matrix\" class")
  
  nrow <- nrow(grid)
  sel <- sample(1:nrow, n, replace = FALSE)

  return(grid[sel,])
}
  
