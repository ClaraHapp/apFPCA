### SRVF transformation approach ###

#' SRVF transformation of warping functions
#'
#' Calculates the SRVF transformation of a set of warping functions
#' \eqn{\gamma}{gamma}: \deqn{\sqrt{\gamma'},}{sqrt(gamma'),} where
#' \eqn{\gamma'}{gamma'} denotes the first derivative of
#' \eqn{\gamma}{gamma}.
#'
#' @param gamma A \code{funData} object containing the warping functions.
#'
#' @return The SRVF transformation of all functions in \code{gamma}, again
#'   as a \code{funData} object.
srvf <- function(gamma)
{
  return(sqrt(diff.funData(gamma)))
}

#' Tangent space transformation of SRVFs
#'
#' This function calculates shooting vectors for SRVFs with respect to a
#' given mean function.
#'
#' @param srvf A \code{funData} object containing the SRVFs.
#' @param mu A \code{funData} object containing the mean function for the
#'   tangent space. Defaults to the constant function, which is associated
#'   with identity warping.
#'
#' @return A \code{funData} object containing the shooting vectors.
tangent <- function(srvf, mu = funData(argvals(srvf), matrix(rep(1, nObsPoints(srvf)), nrow = 1)))
{
  if(nObs(mu) != 1)
    stop("The mean function mu must contain a single function.")
  
  eta <- norm(mu, squared = TRUE)
  
  tmp <- scalarProduct(srvf, mu) / eta
  
  # tmp is guaranteed to lie between 0 and 1, all deviations have numerical reasons
  tmp[tmp < 0] <- 0
  tmp[tmp >= 1] <- 1 - sqrt(.Machine$double.eps)
  
  # angles
  theta <- acos(tmp)
  
  v <- srvf * 0 # initialize return object
  for(i in 1:nObs(srvf))
    v@X[i,] <- theta[i]/(sqrt(eta) * sin(theta[i]))* srvf@X[i,] - cos(theta[i]) * mu@X[1,]
  
  return(v)
}

#' Tangent space transformation for warping functions
#'
#' Combine SRVF transformation with projection on a tangent space for
#' warping functions.
#'
#' @param gamma A \code{funData} object containing the warping functions.
#' @param mu A \code{funData} object containing the mean function for the
#'   tangent space in the warping space. Defaults to identity warping.
#'
#' @return A \code{funData} object containing the shooting vectors.
tangent.warp <- function(gamma, mu = funData(argvals(gamma), matrix(argvals(gamma)[[1]], nrow = 1)))
{
  return(tangent(srvf(gamma) , srvf(mu)))
}

#' Map tangent space vectors to SRVFs
#'
#' This function maps shooting vectors in a tangent space to their SRVFs.
#'
#' @param v A \code{funData} object containing the shooting vectors.
#' @param mu A \code{funData} object containing the mean function for the
#'   tangent space. Defaults to the constant function, which is associated
#'   with identity warping.
#'
#' @return A \code{funData} object containing the srvfs.
tangentInv <- function(v, mu = funData(argvals(v), matrix(rep(1, nObsPoints(v)), nrow = 1)))
{
  if(nObs(mu) != 1)
    stop("The mean function mu must contain a single function.")
  
  eta <- norm(mu)
  
  srvf <- v*0 # initialize result object
  nv <- norm(v, squared = FALSE)
  
  for(i in 1:nObs(v)) 
    srvf@X[i,]   <- sqrt(eta) * sin(nv[i])/nv[i]*v@X[i,] + cos(nv[i]) * mu@X[1,]
  
  return(srvf)
}

#' Inverse SRVF transformation to warping functions
#'
#' Calculates the inverse SRVF transformation of a set of SRVFs.
#' 
#' @param srvf A \code{funData} object containing the SRVFs.
#'
#' @return The inverse SRVF transformation of all functions in \code{srvf}, again
#'   as a \code{funData} object.
srvfInv <- function(srvf){
  return(cumInt.funData(srvf^2))
}

#' Inverse tangent space transformation for warping functions
#'
#' Combine projection from a tangent space to SRVF and inverse SRVF
#' transformation to warping functions.
#'
#' @param v A \code{funData} object containing the functions in tangent space.
#' @param mu A \code{funData} object containing the mean function for the
#'   tangent space in the warping space. Defaults to identity warping.
#'
#' @return A \code{funData} object containing the warping functions.
tangentInv.warp <- function(v, mu = funData(argvals(v), matrix(argvals(v)[[1]], nrow = 1)))
{
  return(srvfInv(tangentInv(v , srvf(mu))))
}


### clr transformation approach ###

#' Centred log-ratio transform for warping functions
#' 
#' @section Warning: The function does not check if the argument is really a warping function!
#' 
#' @param f A funData object that represents a warping function. 
#' 
#' @result A funData object, which is the cenred log-ratio transform (CLR) applied to the derivative of f.
clr.warp <- function(f)
{
  clr(diff.funData(f))
}


#' Centred log-ratio transform
#' 
#' Centred-log ratio transform for functional data objects. The domain needs not necessarily be [0,1].
#' 
#' @section Warning: it is not checked if g is positive and integrable!
#' 
#' @param g A funData object
#' 
#' @return The CLR of g
clr <- function(g)
{
  if(any(g@X < 0, na.rm = TRUE))
  {
    warning("Negative values found, set to NA")
    g@X[g@X < 0] <- NA
   intG <- integrate(as.irregFunData(log(g))) # do integration on irregFunData object
  }
  else
  {
    intG <- integrate(log(g))
  }
    
  return(log(g) - 1/diff(range(g@argvals[[1]])) * intG)
}

#' Inverse centred log-ratio transform
#' 
#' @param g A funData object
#' 
#' @return A funData object, that corresponds to the inverse CLR of g.
clrInv <- function(g)
{
  # map from L^2 to B^2 (densities)
  return(exp(g)/integrate(exp(g)))
}

#' Inverse CLR for warping functions
#'
#' The function calculates an inverse CLR for a funData object g (which gives
#' a density) and integrates it to a warping function. In particular, if the
#' density maps from [a,b], the resulting warping function h has h(a) = a and
#' h(b) = b.
#' 
#' @param g A funData object
#' 
#' @return A funData object corresponding to the associated warping function.
clrInv.warp <- function(g)
{
  h <- clrInv(g)
  # integrate densities to obtain warping functions
  # warping functions h fulfill h(a) = a, h(b) = b
  return(min(h@argvals[[1]]) + diff(range(h@argvals[[1]]))  * cumInt.funData(h))
}


### Petersen and Mueller (2016) transformations ###

#' Log-hazard-transformation
#'
#' Log-hazard-transformation for density functions as defined in Petersen
#' & Mueller (2016), including differentiation of warping functions.
#'
#' @param g A \code{funData} object containing the warping functions.
#' @param delta Cutoff for calculating the hazard function. 
#'   Defaults to \code{0.05}.
#'
#' @return A \code{funData} object containing the transformed functions.
LH.warp <- function(g, delta = 0.05)
{
  # extract argvals
  argvals <- g@argvals[[1]]
  # renormalize to cdf, mapping to [0,1]:
  gN <- (g - min(argvals)) / diff(range(argvals))
  
  # calculate density
  f <- diff.funData(gN) 
  
  gN@X[, argvals > max(argvals) - delta * diff(range(argvals))] <- NA
  
  # trafo
  v <- log(f / (1 - gN))

  return(v)
}  


#' Inverse Log-hazard-transformation
#'
#' Inverse log-hazard-transformation for density functions as defined in Petersen
#' & Mueller (2016), including differentiation of warping functions.
#'
#' @param v A \code{funData} object containing the functions in L2.
#' @param delta Cutoff for calculating the hazard function.
#'   Defaults to \code{0.05}.
#'   
#' @return A \code{funData} object containing the warping functions.
LHInv.warp <- function(v, delta = 0.05)
{
  argvals <- v@argvals[[1]]
  thresh <- max(argvals) - delta * diff(range(argvals))
  v1 <- extractObs(v, argvals = argvals[argvals < thresh])
  
  # density level
  f1 <- exp(v1 - cumInt.funData(exp(v1)))
  f2 <- 1/delta * exp(-1*integrate(exp(v1)))
  
  f <- funData(argvals = argvals, cbind(f1@X, t(sapply(f2, rep, each = length(which(argvals >= thresh))))))
  
  # warping level
  g <- min(argvals) + cumInt.funData(f) * diff(range(argvals))
  
  return(g)
}  
  
#' Log-quantile-transformation
#'
#' Log-quantile-transformation for density functions as defined in
#' Petersen & Mueller (2016), including differentiation of warping
#' functions.
#'
#' @param g A \code{funData} object containing the warping functions.
#' 
#' @return A \code{funData} object containing the transformed functions.
LQ.warp <- function(g)
{
  # quantile function
  Q = invert.funData(g)
  
  # renormalize to cdf, mapping to [0,1]:
  gN <- (g - min(g@argvals[[1]])) / diff(range(g@argvals))
  
  # calculate density
  f <- diff.funData(gN) 
  
  # trafo
  v <- -1 * log(warp.funData(f, Q, smooth = TRUE))
  
  return(v)
}


#' Inverse Log-quantile-transformation
#'
#' Inverse log-quantile-transformation for density functions as defined in
#' Petersen & Mueller (2016), including retransformation to warping
#' functions.
#'
#' @param v A \code{funData} object containing the functions in L2.
#' 
#' @return A \code{funData} object containing the warping functions.
LQInv.warp <- function(v)
{
  # Quantile function
  Q <- cumInt.funData(exp(v)) / integrate(exp(v))
 
  # Warpign function is the inverse of Q, appropriately scaled
  g <- min(v@argvals[[1]]) + invert.funData(Q) * diff(range(v@argvals))
  
  return(g)
}



### Others ###

#' Discrete differentiation for funData objects
#'
#' For the inner points (2:(nObsPoints - 1)), central differentiation is
#' used, which borrows information from the preceding (i-1) as well as the
#' following (i+1) value for calculating the gradient. For the boundary
#' values (1, nObsPoints(f)) only one neighbouring value is used: Forward
#' differentiation for the left bound, backward differentiation for the right
#' bound.
#'
#' @param f A funData object. Must have a one-dimensional domain, otherwise an error is thrown.
#'
#' @return Another funData object, which corresponds to the derivative of f.
diff.funData <- function(f)
{
  if(dimSupp(f) > 1)
    stop("Implementation is only for one-dimensional domains.")
  
  nP <- nObsPoints(f)
  x <- argvals(f)[[1]]
  
  g <- array(NA, dim = dim(f@X)) # initialize g matrix for differentiation
  
  # left bound: forward differentiation
  g[,1] <- (f@X[,2] - f@X[,1]) / truncX(x[2] - x[1])
  
  # right bound: backward differentiation
  g[,nP] <- (f@X[,nP] - f@X[,nP - 1]) / truncX(x[nP] - x[nP - 1])
  
  # all other points: central differentiation
  g[,2:(nP-1)] <- (f@X[,3:nP] - f@X[,1:(nP-2)]) / truncX(x[3:nP] - x[1:(nP-2)])
  
  return(funData(x, g))
}

#' Cumulative integration of funData objects
#'
#' Calculate the integral over a funData object from zero to all points in
#' the observation grid. For numerical stability, the calculation is made
#' backward for lower values and forward for higher values. This avoids
#' integrating over only a few points.
#'
#' @param h A funData object. Must have a one-dimensional domain, otherwise
#'   an error is thrown.
#'
#' @return A funData object, containing the integrated values for each
#'   observation point.
cumInt.funData <- function(h)
{
  if(dimSupp(h) > 1)
    stop("Implementation is only for one-dimensional domains.")
  
  x <- argvals(h)[[1]]
  intH <- integrate(h) # integral over the full domain, needed for calculating the integrals "backward" for lower values
  
  # split integration into two parts to guarantee stability at the boundaries
  n <- max(which(x < median(x)))
  
  cumInt <- lapply(1:nObsPoints(h), function(ind){
    if(ind <= n)
      intH - integrate(extractObs(h, argvals = x[ind:nObsPoints(h)]))
    else
      integrate(extractObs(h, argvals = x[1:ind]))
  })
  
  # cbind & list can handle all numbers of observations (incl. 1)
  return(funData(h@argvals, do.call("cbind", cumInt)))
}


#' Inverse of funData
#'
#' This function returns the inverse of a funData object
#' representing a set of functions on a common grid. The result is again a funData
#' object. The method is based on the smooth.spline function in R (stats)
#' 
#' @param f The functions to be inverted, passed as a funData object
#' @param ... Options to be passed to smooth.spline
invert.funData <- function(f,...)
{
  argvals <- f@argvals[[1]]
  
  return(funData(argvals = argvals, 
          X = t(apply(f@X, 1, function(z){mgcv::predict.gam(mgcv::gam(argvals ~ s(z)), newdata = data.frame(z = argvals))}))))
}


#' Inverse of warping function
#'
#' This function returns the inverse of a funData object representing a
#' set of warping functions on a common grid. The result is again a
#' funData object and respects the monotonicity of the input functions.
#' The method is based on the splinefun function in R (stats)
#'
#' @param f The functions to be inverted, passed as a funData object
#' @param ... Options to be passed to splinefun
invert.warps <- function(f,...)
{
  argvals <- f@argvals[[1]]
  
  return(funData(argvals = argvals, 
                 X = t(apply(f@X, 1, function(z){
                   splinefun(x = z, y = argvals, method = "monoH.FC")(argvals)
                 }))))
}


#' Discrete warping
#' 
#' The calculate the warping of a funData object represented by a funData object containing warping functions.
#' 
#' @section Warning: The function does not check if all elements of w are correct warping functions.
#' 
#' @param f The funData object to be warped
#' @param w The funData object containing the warping functions.
#' @param smooth Logical, should the warped functions be smoothed (via gam?)
#' 
#' @return A funData object containing the warped functions.
warp.funData <- function(f, w, smooth = FALSE)
{
  # check if functions have the same domain
  if(! all.equal(range(argvals(f)), range(argvals(w))))
    stop("f and w need to have the same domain!")
  
  res <- matrix(NA, nrow = nObs(f), ncol = length(w@argvals[[1]])) # initialize result (unwarped functions)
  
  # for all functions in f
  for(i in 1:nObs(f))
  {
    allDiff <- outer(w@X[i,], f@argvals[[1]], function(x,y){abs(x-y)}) # calculate differences between all x values
    p <- f@X[i,apply(allDiff, 1, which.min)] # choose nearest x- value
    
    if(smooth)
    {
      xi <- argvals(w)[[1]]
      GAM <- mgcv::gam(p ~ s(xi, bs = "ps")) 
      p <- predict(GAM)
    }
    
    res[i,] <- p
  }
  
  return(funData(w@argvals,res))
}

#' Truncate very small absolute values to 0
#' 
#' @param x The value to truncate
#' @param eps The threshold for truncating
#' 
#' @return The truncated value (either the original or 0, if it is too small)
truncX <- function(x, eps = sqrt(.Machine$double.eps))
{
  return(ifelse(abs(x) < eps, eps, x))
}

#' Non-smooth FPCA for functions on a regular grid
#'
#' This function calculates a functional PCA for funData objects on a
#' regular grid, without smoothing the mean oder covariance functions.
#'
#' @param f The funData object containg the functions to analyze.
#' @param pve The proportion of variability explained. Defaults to 0.99.
#'
#' @return A list containing the mean function, the principal components,
#'   eigenvalues and scores.
nonSmoothFPCA <- function(f, pve = 0.99)
{
  # Demean
  mu <- meanFunction(f)
  f <- f - mu
  
  # Calculate PCA (analogous to MFPCA::PACE, but no smoothing)
  w <- funData:::.intWeights(argvals(f)[[1]], method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% cov(X(f)) %*% Wsqrt
  evalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues <- replace(evalues, which(evalues <= 0), 0)
  npc <- min(which(cumsum(evalues)/sum(evalues) > pve))
  efunctions <- funData(argvals(f), 
                        t(matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], 
                                 nrow = nObsPoints(f), ncol = npc)))
  evalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  
  scores <- matrix(NA, nrow = nObs(f), ncol = nObs(efunctions))
  for(i in 1:nObs(efunctions))
    scores[,i] <- scalarProduct(f - meanFunction(f), efunctions[i])
  
  return(list(mu = mu,
              functions = efunctions,
              values = evalues,
              scores = scores))
}
