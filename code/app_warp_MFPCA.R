### Calculate SRVF warping & MFPCA ###

### Log-transform smoothed velocities
groundVel_smooth <- funData(argvals = x, 
                       X = as.matrix(log1p(seissol$groundVel_tw_smooth[ind, ])))

### Calculate SRVF warping
SRVFwarp <- fdasrvf::time_warping(f = t(groundVel_smooth@X), time = groundVel_smooth@argvals[[1]], 
                                  lambda = 0,# controls elasticity
                                  method = "mean", # Karcher mean
                                  showplot = FALSE, #do not show plots of functions 
                                  smooth_data = FALSE, #no box filter
                                  MaxItr = 500) # maximum number of iterations

# extract warping results in form of funData objects
h <- funData(groundVel_smooth@argvals, min(groundVel_smooth@argvals[[1]]) + t(SRVFwarp$gam)  * diff(range(groundVel_smooth@argvals)))
gamma <- invert.warps(h) # invert, for our notation
aligned <- funData(groundVel_smooth@argvals, t(SRVFwarp$fn))

### Calculate MFPCA

# Create multiFunData object
m <- multiFunData(clr.warp(gamma), aligned)

# univariate FPCA (use as given basis for faster calculation)
pca <- list(MFPCA::PACE(m[[1]]),
            nonSmoothFPCA(m[[2]]))
uniEx <- list(list(type = "given", functions = pca[[1]]$functions, scores = pca[[1]]$scores),
              list(type = "given", functions = pca[[2]]$functions, scores = pca[[2]]$scores))

# subtract univariate means
mu <- multiFunData(pca[[1]]$mu, pca[[2]]$mu)
m <- m - mu

# Optimize weight of warping element
findWeight <- function(C, M){
  if(options()$verbose)
    cat("C: ", C, "\n")
  
  PCAm <- MFPCA::MFPCA(m, M = M, 
                       uniExpansions = uniEx,
                       weights = c(C,1), fit = TRUE)
  
  # reconstruction
  xHat <- warp.funData(mu[[2]] + PCAm$fit[[2]], clrInv.warp(mu[[1]] + sqrt(C)*PCAm$fit[[1]]), smooth = FALSE) # really smooth?
  
  mean(norm(xHat - groundVel_smooth))
}

bestApprox <- optimize(findWeight, interval = c(0,10), M = 10)

# save all
save(seissol, groundVel_smooth, SRVFwarp, h, aligned, gamma, mu, m, bestApprox, file = "../data/SRVF_seis_tw_smooth_rawPCA.Rdata")