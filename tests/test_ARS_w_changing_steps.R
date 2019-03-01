#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Run this to test how changes in scale steps affect the FPT output ###
## Data here have been run through FULL_ANALYSIS_TEST until the scaleARS() step
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

steps <- c(0.5, 1, 5, 10, 20, 40) ## set scale steps to test (km)

scales.list <- list()
scales.list[[1]] <- seq(0.1, 250, steps[1])
scales.list[[2]] <- seq(0.1, 250, steps[2])
scales.list[[3]] <- seq(0.1, 250, steps[3])
scales.list[[4]] <- seq(0.1, 250, steps[4])
scales.list[[5]] <- seq(0.1, 250, steps[5])
scales.list[[6]] <- seq(0.1, 250, steps[6])
# scales.list[[7]] <- c(seq(10, 25, 1), seq(30, 50, 5), 75, seq(100, 250, 50))


ARS.list <- list()
# path <- 'figures/scaleARS'  # set plot path and parameters # set path if you want to save all FPT output figs. for comparison

for(i in 1:length(scales.list)){
  
  Scales <- scales.list[[i]]
  
  # png(file.path(path, paste(steps[i], "km", ".png", sep="")), width=600, height = 500, res=100) 
  fpt.scales <- scaleARS(DataGroup=Trips, Scales=Scales, Peak="Flexible") 
  # dev.off()
  
  ARS.list[[i]] <- fpt.scales
  
}

ARS.list

ARS.df <- data.frame(ARS=do.call("rbind", ARS.list))
ARS.df$step <- steps

# png(file.path(path, paste("WAAL_1-40km_steps_ARS", ".png", sep="")), width=600, height = 500, res=100)
plot(ARS.df$step, ARS.df$ARS, type='b', xlab="ARS step sequence (km)", ylab="median ARS from FPT (km)", ylim=c(10, 250))
# dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Vary the starting position of the scale range #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ARSn.list <- list()

for(s in 0:9) { ## e.g. adding 0-9km to the starting scale

scales.list2 <- lapply(scales.list, function(x) x+s)

  for(i in 1:length(scales.list2)){
    
    Scales <- scales.list2[[i]]
    
    # png(file.path(path, paste(steps[i], "km", ".png", sep="")), width=600, height = 500, res=100)
    fpt.scales <- scaleARS(DataGroup=DataGroupTrips, Scales=Scales, Peak="Flexible") 
    # dev.off()
    
    ARS.list[[i]] <- fpt.scales
    
  }

ARS.list

ARS.df <- data.frame(ARS = do.call("rbind", ARS.list))
ARS.df$step <- steps
ARS.df$start_scale <- rep(10 + s)
ARS.df$ARS <- ARS.df$ARS + (10 - ARS.df$start_scale) ## Adjust the final ARS scale by the 
# png(file.path(path, paste("WAAL_1-40km_steps_ARS", ".png", sep="")), width=600, height = 500, res=100)
# plot(ARS.df$step, ARS.df$ARS, type='b', xlab="ARS step sequence (km)", ylab="median ARS from FPT (km)", ylim=c(10, 250))
# dev.off()

ARSn.list[[s+1]] <- ARS.df

}

path <- 'figures/scaleARS/shift_range'  # set plot path and parameters

png(file.path(path, paste("shift start_by10", ".png", sep="")), width=600, height = 500, res=100)
with(ARSn.list[[1]], plot(step, ARS, type='b', xlab="ARS step sequence (km)", ylab="median ARS from FPT (km)", ylim=c(10, max(ARS))))
with(ARSn.list[[2]], points(step, ARS, type='b', col=2))
with(ARSn.list[[3]], points(step, ARS, type='b', col=3))
with(ARSn.list[[4]], points(step, ARS, type='b', col=4))
with(ARSn.list[[5]], points(step, ARS, type='b', col=5))
with(ARSn.list[[6]], points(step, ARS, type='b', col=6))
with(ARSn.list[[7]], points(step, ARS, type='b', col=7))
with(ARSn.list[[8]], points(step, ARS, type='b', col=8))
with(ARSn.list[[9]], points(step, ARS, type='b', col=9))
with(ARSn.list[[10]], points(step, ARS, type='b', col=10))
with(ARSn.list[[1]], points(step, ARS, type='b', col=1, cex=1.5, pch=20))
text(x=20, y=max(ARSn.list[[1]]$ARS), "10:20 - 250:260, by=c(1, 2, 5, 10, 20, 40km)")
dev.off()
