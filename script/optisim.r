# This file is part of the inspectr package. The license of that package applies.

if (!isGeneric("optisim"))
  setGeneric("optisim", function(x, size, ...)
    standardGeneric("optisim"))

setMethod("optisim", "Spectra", 
  function(x, size, ...) {
    idx <- optisim.matrix(spectra(x), size = size, ...)
    x[idx,]
  }
)


# Distance of each point of a matrix x to
# a unique value
.distToMean <- function(x, mn){
  apply(x, 1, function(x) dist(rbind(mn, x)))
}

# Distance of each point in a matrix x to
# each point in another matrix sub
.distToSubset <- function(x, sub) {
  apply(x, 1, function(y) apply(sub, 1, function(x) dist(rbind(x, y))))
}

## OptiSim, an efficient implementation of the KS algorithm
##
## -- use it on PCs to compress the data !! --
##
#
# In my (short) experience, optisim() is usually ~10 times
# faster than kenstone().
#
# B is the number of objects in the random subset
#
# idx <- optisim(spectra(foo), size = 15)
# idx <- optisim(foo[, 'PC1', 'PC2', 'PC3'], size = 15)
# 
optisim.matrix <- function(x, size, B = round(0.10*nrow(x)), progress = TRUE, ...) { 
  
  n <- nrow(x)
  m <- ncol(x)
  min_x <- apply(x, 2, min, na.rm = TRUE)
  max_x <- apply(x, 2, max, na.rm = TRUE)
  range <- max_x - min_x
  V <- cumprod(range)[length(range)]
  Vs <- (1/size)*V
  # to avoid warning on +Inf values of Gamma
  ifelse(m > 100, G <- Inf, G <- gamma(m/2 + 1)) 
  R <- (Vs/(sqrt(pi^m)/G))^(1/m)

  mean_x <- colMeans(x)
  d <- colSums((x - matrix(rep(mean_x, n), nrow = n, byrow = TRUE))^2)
  i <- min(d)
  r <- which.min(d)
  
  model <- r
  A <- 1:n
  A <- A[-r]
  recycling <- NULL
  mS <- 1
  ma <- length(A)

  if (progress) 
    pb <- txtProgressBar(min = 1, max = size, style = 3)

  while(mS < size) {

    r <- sample(1:ma, size = ma)

    if (B < ma) 
      r <- r[1:B]
    
    odleglosc <- .distToSubset(x[model, , drop = FALSE], x[A[r], , drop = FALSE])
    
    if (mS > 1) 
      odleglosc <- min(odleglosc) # apply(odleglosc, 2, min)

    remove <- which(odleglosc < R)
    i1 <- max(odleglosc)
    i2 <- which.max(odleglosc)
    model <- c(model, A[r[i2]])

    if (length(remove) > 0) 
      A <- A[-r[remove]]
    
    A = A[-r[i2]] # GEm: make sure we remove used samples from the population
    r <- r[-i2]
    recycling <- c(recycling, r)

    if (length(A) == 0)
      A <- unique(recycling)
    
    ma <- length(A)
    mS <- length(model)

    if (progress)
      setTxtProgressBar(pb, mS)
  }

  if (progress)
    close(pb)

  model
}
