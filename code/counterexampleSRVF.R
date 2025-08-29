### Counterexample
library("funData")
source("utils.R")
set.seed(4321)

# Simulate warping functions
N <- 50
ki <- rgamma(N,5,5)

t <- seq(0.001,1-0.001, 0.005)
gamma <- funData(t, t(sapply(ki, function(j){t^j})))
gamma2 <- funData(t, matrix(t^{5}, nrow = 1)) # new function for prediction

C <- 1 # weighting factor for PC visualization

### SRVF trafo (tangent space)
x <- tangent.warp(gamma) # transform to tangent space, use default mu

# calculate FPCA
pcaSRVF <- MFPCA::PACE(x)
pcaSRVF$values/sum(pcaSRVF$values) # first PC explains most of variability in L2

# visualize effect of first PC in tangent space
pcsTangent <- pcaSRVF$mu + C*sqrt(pcaSRVF$values[1]) * 
  funData(pcaSRVF$functions@argvals, rbind(pcaSRVF$functions@X[1,], -1 * pcaSRVF$functions@X[1,]))

# sanity checks
integrate(x) # tangent space functions orthogonal to constant functions
norm(srvfInv(pcsTangent)) # on SRVFs of PCs on unit sphere?

# prediction of new curve based on first PC
predSRVF <- tangentInv.warp(pcaSRVF$mu + MFPCA:::expandBasisFunction(scores = matrix(scalarProduct(tangent.warp(gamma2)  - pcaSRVF$mu, pcaSRVF$functions[1]), nrow = 1),
                                                                     functions = pcaSRVF$functions[1]))


##### CLR trafo (Bayes space)
bayes <- clr.warp(gamma) # transform to Bayes space

# calculate FPCA
pcaBayes <- MFPCA::PACE(bayes) 
pcaBayes$values / sum(pcaBayes$values)

# visualize effect of first PC in Bayes space
pcsBayes <- pcaBayes$mu + C*sqrt(pcaBayes$values[1]) *funData(pcaBayes$functions@argvals, rbind(pcaBayes$functions@X[1,], -1 * pcaBayes$functions@X[1,]))

# prediction of new curve based on first PC
predBayes <- clrInv.warp(pcaBayes$mu + MFPCA:::expandBasisFunction(scores = matrix(scalarProduct(clr.warp(gamma2)  - pcaBayes$mu, pcaBayes$functions[1]), nrow = 1),
                                                                   functions = pcaBayes$functions[1]))

### Log-hazard transformation
logHazard <- LH.warp(gamma) # transform to L2 space

# calculate FPCA
pcaLH <- MFPCA::PACE(logHazard) 
pcaLH$values / sum(pcaLH$values)

# visualize effect of first PC in L2 space
pcsLH <- pcaLH$mu + C*sqrt(pcaLH$values[1]) *funData(pcaLH$functions@argvals, rbind(pcaLH$functions@X[1,], -1 * pcaLH$functions@X[1,]))

# prediction of new curve based on first PC
predLH <- LHInv.warp(pcaLH$mu + MFPCA:::expandBasisFunction(scores = matrix(scalarProduct(as.irregFunData(LH.warp(gamma2)  - pcaLH$mu), pcaLH$functions[1]), nrow = 1),
                                                            functions = pcaLH$functions[1]))

### Log-quantile transformation
logQuantile <- LQ.warp(gamma) # transform to L2 space

# calculate FPCA
pcaLQ <- MFPCA::PACE(logQuantile) 
pcaLQ$values / sum(pcaLQ$values)

# visualize effect of first PC in L2 space
pcsLQ <- pcaLQ$mu + C*sqrt(pcaLQ$values[1]) *funData(pcaLQ$functions@argvals, rbind(pcaLQ$functions@X[1,], -1 * pcaLQ$functions@X[1,]))

# prediction of new curve based on first PC
predLQ <- LQInv.warp(pcaLQ$mu + MFPCA:::expandBasisFunction(scores = matrix(scalarProduct(as.irregFunData(LQ.warp(gamma2)  - pcaLQ$mu), pcaLQ$functions[1]), nrow = 1),
                                                            functions = pcaLQ$functions[1]))


#### Plot the results
library(ggplot2)

pdf("counterexample_Plots.pdf")

# SRVF
autoplot(x, col = "grey", alpha = 0.4) +
  autolayer(pcaSRVF$functions, obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(pcsTangent, col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(pcsTangent[1,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(pcsTangent[2,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(tangentInv(x), col = "grey", alpha = 0.4) +
  autolayer(tangentInv(pcaSRVF$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(tangentInv(pcsTangent), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(tangentInv(pcsTangent[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(tangentInv(pcsTangent[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(gamma, col = "grey", alpha = 0.4) +
  autolayer(tangentInv.warp(pcaSRVF$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(tangentInv.warp(pcsTangent),  col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(tangentInv.warp(pcsTangent[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(tangentInv.warp(pcsTangent[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  autolayer(gamma2, col = "darkgreen", lwd = 1, lineend  = "round") + 
  autolayer(predSRVF, col = "darkgreen", lwd = 1, lineend  = "round", lty = 3) + 
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

# Bayes trafo
autoplot(bayes, col = "grey", alpha = 0.4) +
  autolayer(pcaBayes$functions, obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(pcsBayes,  col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(pcsBayes[1,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(pcsBayes[2,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(clrInv(bayes), col = "grey", alpha = 0.4) +
  autolayer(clrInv(pcaBayes$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(clrInv(pcsBayes), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(clrInv(pcsBayes[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(clrInv(pcsBayes[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(gamma, col = "grey", alpha = 0.4) +
  autolayer(clrInv.warp(pcaBayes$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(clrInv.warp(pcsBayes), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(clrInv.warp(pcsBayes[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(clrInv.warp(pcsBayes[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  autolayer(gamma2, col = "darkgreen", lwd = 1, lineend  = "round") + 
  autolayer(predBayes, col = "darkgreen", lwd = 1, lineend  = "round", lty = 3) + 
  labs(x = "", title = "")+
  theme_bw(base_size = 20)


# log-hazard trafo
autoplot(logHazard, col = "grey", alpha = 0.4) +
  autolayer(pcaLH$functions, obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(pcsLH, col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(pcsLH[1,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(pcsLH[2,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  ggplot2::geom_vline(xintercept = 1- 0.05, lwd = 0.5, lty = 3) +
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

# here, we get the warping functions directly and diff them to obtain densities
autoplot(diff.funData(gamma), col = "grey", alpha = 0.4) +
  autolayer(diff.funData(LHInv.warp(pcaLH$functions)), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(diff.funData(LHInv.warp(pcsLH)), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(diff.funData(LHInv.warp(pcsLH[1,argvals = t[seq(1, 200, by = 15)]])), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(diff.funData(LHInv.warp(pcsLH[2,argvals = t[seq(1, 200, by = 15)]])), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  ggplot2::geom_vline(xintercept = 1- 0.05, lwd = 0.5, lty = 3) +
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(gamma, col = "grey", alpha = 0.4) +
  autolayer(LHInv.warp(pcaLH$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(LHInv.warp(pcsLH), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(LHInv.warp(pcsLH[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(LHInv.warp(pcsLH[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  autolayer(gamma2, col = "darkgreen", lwd = 1, lineend  = "round") + 
  autolayer(predLH, col = "darkgreen", lwd = 1, lineend  = "round", lty = 3) + 
  ggplot2::geom_vline(xintercept = 1- 0.05, lwd = 0.5, lty = 3) +
  labs(x = "", title = "")+
  theme_bw(base_size = 20)


# log-quantile trafo
autoplot(logQuantile, col = "grey", alpha = 0.4) +
  autolayer(pcaLQ$functions, obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(pcsLQ, col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(pcsLQ[1,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(pcsLQ[2,argvals = t[seq(1, 200, by = 15)]], geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

# here, we get the warping functions directly and diff them to obtain densities
autoplot(diff.funData(gamma), col = "grey", alpha = 0.4) +
  autolayer(diff.funData(LQInv.warp(pcaLQ$functions)), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(diff.funData(LQInv.warp(pcsLQ)), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(diff.funData(LQInv.warp(pcsLQ[1,argvals = t[seq(1, 200, by = 15)]])), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(diff.funData(LQInv.warp(pcsLQ[2,argvals = t[seq(1, 200, by = 15)]])), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  labs(x = "", title = "")+
  theme_bw(base_size = 20)

autoplot(gamma, col = "grey", alpha = 0.4) +
  autolayer(LQInv.warp(pcaLQ$functions), obs = 1, col = "black", lwd = 1.5, lineend  = "round")+
  autolayer(LQInv.warp(pcsLQ), col = "darkblue", size = 0.5, alpha = 0.5, lineend  = "round")+
  autolayer(LQInv.warp(pcsLQ[1,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "+", size = 5, stroke = 3)+
  autolayer(LQInv.warp(pcsLQ[2,argvals = t[seq(1, 200, by = 15)]]), geom = "point",  col = 4, shape = "-", size = 5, stroke = 3)+
  autolayer(gamma2, col = "darkgreen", lwd = 1, lineend  = "round") + 
  autolayer(predLQ, col = "darkgreen", lwd = 1, lineend  = "round", lty = 3) + 
  labs(x = "", title = "")+
  theme_bw(base_size = 20)
dev.off()
