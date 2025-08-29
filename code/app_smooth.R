# Parallelize smoothing
plan(multiprocess(workers = 30))
seissol = readRDS("../data/data_northridge.rds")

# Preprocess
data = seissol$Bodenbewegung %>% 
  mutate(id = seissol$Seismogram:seissol$Simulation) %>% 
  gather(key = "time", value = "y", -id) %>% 
  mutate(time = as.numeric(str_extract(time, "[0-9\\.]+$")), 
         y = zapsmall(y)) %>% 
  nest(-id, .key = "boden")

# Function to smooth one curve using the Tweedie-distribution
# we have true zeroes so use tw, p < 2 this which has point mass at 0. tw(link = "log")
smooth_one = function(d, k = 40, bs = "cr", m = 2, family = tw(), optimizer = c("outer","newton"), 
                      sp = NULL, return = c("fitted.values", "sp"), offset=NULL) {
  m = gam(y ~ s(time, k = k, bs = bs, m = m), offset = offset,
          family = family, optimizer = optimizer, sp = sp, data = d)
  return(m[return])
}

# Smooth ground velocity using smooth_one function
groundVel_tw_smooth = 
  future_map(data$boden,
             ~ cbind(.x, fit = smooth_one(.x, sp = 2.7)[["fitted.values"]]))

# Add id
data_groundVel_tw_smooth = groundVel_tw_smooth %>% 
  imap(~cbind(id = .y, .x[, c(1,3)])) %>% #drop original data, add id
  bind_rows %>% 
  spread(key = time, value = fit)

# Add smoothed data to data frame
seissol$groundVel_tw_smooth = data_groundVel_tw_smooth[,-1]

# Save data frame
saveRDS(seissol, "../data/seissol_tweedie_smooth.rds")