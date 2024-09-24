library(tidyverse)
library(sizeSpectra)

source("landreth/bin_data.R")
ldat <- readRDS("landreth/landreth_fishmacros_data.rds")
ldat

llist <- group_split(ldat, stream)
length(llist)
stream_names <- unique(ldat$stream)


peak_edge <- function(dat_list){
  
  llist <- dat_list
  
  out <- data.frame(
    stream = vector("character", length(llist)),
    count_min = vector("numeric", length(llist)),
    sum_min = vector("numeric", length(llist)))
  out
  
  for(i in 1:length(llist)){
    # read in one data set
    dat <- llist[[i]]
    # save stream name in output
    out$stream[i] <- unique(dat$stream)
    dat <- dat |>
      select(x = dw_g, 
             counts)
    x_obs_bin <- bin_data(counts = dat,
                          binWidth = "2k")
    # calculate minimum values
    sum_peak <- which.max(x_obs_bin$binVals$binSumNorm)
    count_peak <- which.max(x_obs_bin$binVals$binCountNorm)
    print(c(sum_peak, count_peak))
    
    # keeps the maximum edge of the bin peak
    # conservative estimate of where the peak must start
    out$sum_min[i] <- x_obs_bin$binVals$binMax[sum_peak]
    out$count_min[i] <- x_obs_bin$binVals$binMax[count_peak]
  }
  return(out)
}

fish_macros_out <- peak_edge(llist)

saveRDS(fish_macros_out, "landreth/peak_min_sizes_fish_macros.RDS")


fish <- readRDS("landreth/landreth_fish_data.rds")
macros <- readRDS("landreth/landreth_macros_data.rds")

fishlist <- group_split(fish, stream)
macroslist <- group_split(macros, stream)
fish_out <- peak_edge(fishlist)
macros_out <- peak_edge(macroslist)

saveRDS(fish_out, "landreth/peak_min_sizes_fish.RDS")
saveRDS(macros_out, "landreth/peak_min_sizes_macros.RDS")
