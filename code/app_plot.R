### Final analysis and plotting ###

### Identify special functions
indMin <- which.min(seissol$hypo.dist[ind])
indMed <- which.min(abs(seissol$hypo.dist[ind] - median(seissol$hypo.dist[ind])))
indMax <- which.max(seissol$hypo.dist[ind])

### Start plotting
pdf(paste("seisPlots_tw_smooth_rawPCA.pdf", sep = ""), height = 5.5, width = 7.5)

### Display original functions and warping results
pOriginal <- autoplot(groundVel_smooth, alpha = 0.05) +
  autolayer(groundVel_smooth, obs = indMin, size = 1.5, col = 2) + # curve with minimal distance
  autolayer(groundVel_smooth, obs = indMed, size = 1.5, col = 3) + # curve with median distance
  autolayer(groundVel_smooth, obs = indMax, size = 1.5, col = 4) + # curve with maximal distance (in this subset)
  autolayer(raw, obs = indMin, alpha = 0.75,  col = 2, geom = "point") + # curve with minimal distance
  autolayer(raw, obs = indMed, alpha = 0.75,  col = 3, geom = "point") + # curve with median distance
  autolayer(raw, obs = indMax, alpha = 0.75,  col = 4, geom = "point") + # curve with maximal distance (in this subset)
  labs(x = "Time", title = "") + 
  theme_bw(base_size = 20)

pAligned <- autoplot(aligned, alpha = 0.05) +
  autolayer(aligned, obs = indMin, size = 1.5, col = 2) + # curve with minimal distance
  autolayer(aligned, obs = indMed, size = 1.5,col = 3) + # curve with median distance
  autolayer(aligned, obs = indMax,size = 1.5, col = 4) + # curve with maximal distance (in this subset)
  labs(x = "Time", title = "")+
  scale_y_continuous(limits = range(groundVel_smooth@X))+
  theme_bw(base_size = 20)

pWarp <- autoplot(gamma, alpha = 0.05) + 
  autolayer(gamma, obs = indMin, size = 1.5, col = 2) + # curve with minimal distance
  autolayer(gamma, obs = indMed, size = 1.5, col = 3) + # curve with median distance
  autolayer(gamma, obs = indMax, size = 1.5, col = 4) + # curve with maximal distance (in this subset)
  labs(x = "Individual Time", y = "Absolute Time", title = "") + 
  theme_bw(base_size = 20)

print(pOriginal)
print(pAligned)
print(pWarp)


### MFPCA based on CLR trafo

# Optimal weight
cat("Optimal weight C: ", round(bestApprox$minimum,2), "\n")

# Rerun FPCA with optimal weight
pca <- list(MFPCA::PACE(m[[1]]),
            MFPCA::PACE(m[[2]]))
uniEx <- list(list(type = "given", functions = pca[[1]]$functions, scores = pca[[1]]$scores),
              list(type = "given", functions = pca[[2]]$functions, scores = pca[[2]]$scores))
PCAopt <- MFPCA::MFPCA(m, M = 10, 
                       uniExpansions = uniEx,
                       weights = c(bestApprox$minimum,1), fit = TRUE)

# reconstruction
xHat <- warp.funData(mu[[2]] + PCAopt$fit[[2]], clrInv.warp(mu[[1]] + sqrt(bestApprox$minimum)*PCAopt$fit[[1]]), smooth = FALSE)
summary(norm(groundVel_smooth - xHat)) # reconstruction accuracy
summary(norm(groundVel_smooth - xHat)/norm(groundVel_smooth))


### Principal component effects

# For first 5 PCs (not all shown in the paper)
for(K in 1:5)
{
  # PC effect
  p.all <- p.phase <- p.amplitude <- autoplot(warp.funData(mu[[2]],  clrInv.warp(mu[[1]]), smooth = FALSE), col = 1, lwd = 1.5) +
    labs(x = "Time", title = paste("PC", K)) + 
    theme_bw(base_size = 20)
  
  for(v in c(-1,1))
  {
    # Calculate PC effect
    PCeffect <- function(alphaWarp, alphaAligned)
    {
      warp.funData(mu[[2]] + v * alphaAligned * sqrt(PCAopt$values[K])*PCAopt$functions[[2]][K], 
                   clrInv.warp(mu[[1]] + v * alphaWarp * sqrt(PCAopt$values[K]) * sqrt(bestApprox$minimum)* PCAopt$functions[[1]][K]), 
                   smooth = FALSE)
    }
    
    p.all <- p.all + 
      autolayer(PCeffect(alphaWarp = 1, alphaAligned = 1), 
                geom = "point",  col = 4, shape = ifelse(v > 0, "+", "-"), size = 5, stroke = 3) +
      autolayer(PCeffect(alphaWarp = 1, alphaAligned = 1),  
                geom = "line",  col = "darkblue", size = 0.5, alpha = 0.5) +
      labs(subtitle = "Full variation",
           title = paste("PC ", K, " (", round(100 * PCAopt$values[K] / sum(PCAopt$values), 2), "% Fréchet variance explained)", sep = ""))
    
    p.phase <- p.phase + 
      autolayer(PCeffect(alphaWarp = 1, alphaAligned = 0), 
                geom = "point",  col = 4, shape = ifelse(v > 0, "+", "-"), size = 5, stroke = 3) +
      autolayer(PCeffect(alphaWarp = 1, alphaAligned = 0),  
                geom = "line",  col = "darkblue", size = 0.5, alpha = 0.5) +
      labs(subtitle = "Phase variation")
    
    p.amplitude <- p.amplitude +
      autolayer(PCeffect(alphaWarp = 0, alphaAligned = 1), 
                geom = "point",  col = 4, shape = ifelse(v > 0, "+", "-"), size = 5, stroke = 3) +
      autolayer(PCeffect(alphaWarp = 0, alphaAligned = 1),  
                geom = "line",  col = "darkblue", size = 0.5, alpha = 0.5) +
      labs(subtitle = "Amplitude variation")
  } 
  
  print(p.all)
  print(p.phase)
  print(p.amplitude)
  
  # spatial distribution
  scorePlot <- ggplot(data.frame(lat = seissol$lat[ind] / 1000, lon = seissol$lon[ind] / 1000, s1 = PCAopt$scores[,K]),
                      aes(x = lon, y = lat, color = s1)) +
    geom_point(size = 3.5) +  scale_color_distiller(palette = 16, name = paste("Scores PC", K)) +
    geom_point(size = 3.6, aes(x = seissol$lon[ind][indMin] / 1000, y = seissol$lat[ind][indMin] / 1000), shape=21, fill = NA, col = "red")+
    labs(x = "Longitude [km]", y = "Latitude [km]", title = "Spatial representation of scores") +
    scale_x_reverse() + # higher values in long = longer distance from Greenwhich = more west
    theme_bw(base_size = 20) 
  
  print(scorePlot)
}

### Score scatterplot (PC 1 & PC 2)
scatter_1_2 <- ggplot(data = data.frame(score1 = PCAopt$scores[,1],
                                        score2 = PCAopt$scores[,2], dist = seissol$hypo.dist[ind] / 1000),
                      aes(x = score1, y = score2, col = dist)) +
  geom_point(alpha = 0.7) +
  scale_color_distiller(palette = 16, name = "Hypocentral\nDistance [km]") +
  labs(x = "Scores PC 1", y = "Scores PC 2", title = "Score scatterplot") +
  theme_bw(base_size = 20)
print(scatter_1_2)

dev.off()
