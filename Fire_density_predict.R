#------------- DESCRIPTION-------------
# The code and data files here accompany the manuscript:
# “An integrated occupancy and home-range model to predict abundance of a wide-ranging, territorial vertebrate” 
# by Morgan W. Tingley, Robert L. Wilkerson, Christine A. Howell, and Rodney B. Siegel.
# Any questions should be directed to Morgan Tingley at <morgan.tingley@uconn.edu>
#  
# This code applies the Black-backed Woodpecker abundance model to a fire with unknown woodpecker density
#
# Output: GIS raster layers of (1) modal BBWO density (pairs/pixel) and (2) standard deviation of density
#--------------------------

## Load required libraries ('rgdal' may require additional installations if computer not previously set-up for R-based GIS)
library(raster)
library(rgdal)
library(R2jags)  # not required but will give errors otherwise

#------------- Import data-------------
## Load Occupancy data and model posterior
folder <- getwd()
setwd("./Model_Occupancy/")
load("Occupancy_data.Rdata")
load("Occupancy_posterior.Rdata") 
# Load Telemetry model posterior 
setwd("../")
setwd("./Model_HomeRange//")
load("HR_posterior.Rdata")
# Load Snag model posterior
setwd("../")
setwd("./Model_Snag//")
load("snag_posterior.Rdata")
load("snag_variable_means.Rdata")

# Load in GIS layers for fire for which you wish to model abundance
# Here, using the example of the Reading fire
# All GIS layers should be appropriately pre-processed: all raster, all same resolution, all same projection
setwd(folder)
setwd("./GIS_Reading")
fire_cc            <- raster("fire_cc.tif")  # Change in canopy cover post-fire (%)
fire_precc         <- raster("fire_precc.tif")  # Pre-fire canopy cover (%)
fire_dem           <- raster("fire_dem_aea.tif")  # Elevation (m)
fire_lat           <- raster("fire_lat.tif")  # Latitude (°)
fire_whr           <- raster("fire_whrtype.tif")  # WHR-forest type (integer links to crosswalk classes)
fire_size          <- raster("fire_whrsize.tif")  # WHR-forest tree size (0-6)
fire_whr_crosswalk <- read.csv("whr_crosswalk.csv")
fire_outline       <- readOGR(dsn = getwd(), layer="fire_outline")  # Fire outline is not necessary but can be used as a boundary for calculating total abundance within fire
setwd(folder)


#------------- Scaling variables-------------
# Most covariates were scaled prior to use in snag and occupancy models. GIS-layers for predicted fire need to be put on the same scale
fire_size_raw     <- fire_size 
elev.res.mod      <- coef(lm(occ_data$elev.raw ~ occ_data$lat.raw))
fire_elev_resid   <- fire_dem - (fire_lat * elev.res.mod[2] + elev.res.mod[1])
fire_elev.st      <- (fire_elev_resid - mean(data.frame(residuals(lm(occ_data$elev.raw ~ occ_data$lat.raw)))[,1])) / 
                     sd(data.frame(residuals(lm(occ_data$elev.raw ~ occ_data$lat.raw)))[, 1])
fire_elev.st2     <- fire_elev.st * fire_elev.st
fire_lat.st       <- (fire_lat - mean(occ_data$lat.raw)) / sd(occ_data$lat.raw)
values(fire_size) <- ifelse(values(fire_size) > 3, 1, 0)
fire_cc.st        <- (fire_cc - mean(occ_data$cc.raw)) / sd(occ_data$cc.raw)
fire_precc.st     <- (fire_precc - mean(occ_data$precc.raw)) / sd(occ_data$precc.raw)
whr_index         <- levels(occ_data$whrtype)
fire_veg          <- fire_whr_crosswalk[values(fire_whr), 3]
fire_veg2         <- match(fire_veg, whr_index)
fire_veg2[is.na(fire_veg2)] <- which(whr_index == "OTH")
# Scaling of values just for snag model
snag_fire_cc.st    <- (fire_cc - snag_data$mu.cc) / snag_data$sd.cc 
snag_fire_cc.st2   <- snag_fire_cc.st * snag_fire_cc.st
snag_fire_precc.st <- (fire_precc - snag_data$mu.precc) / snag_data$sd.precc

# Trying to apply the model to all pixels in a fire at 1 time can be computationally time consuming and space inefficient.
# This code has been written to analyze one pixel at a time, and to loop the model over all pixels.
# For efficiency, not all pixels (including NA's) are included in the loop, only those with WHR values
# Note: for greater efficiency, smaller subsets could be used
fire_outline_rast <- rasterize(x = fire_outline, y = fire_cc)
in_fire           <- (1 - is.na(values(fire_outline_rast)))
is_conifer        <- (1 - is.na(values(fire_whr)))
pix               <- intersect(which(in_fire == 1), which(is_conifer == 1))

#------------- Running model-------------
`expit` <- function(x) {exp(x) / (1 + exp(x))}
Mode    <- function(z) {density(z)$x[which.max(density(z)$y)]}

# Pull out posteriors from which to sample parameter space
post_psi  <- occ_post$BUGSoutput$sims.list
post_hr   <- hr_post$BUGSoutput$sims.list
post_snag <- snag_post$BUGSoutput$sims.list

# Extract values from rasters (for speed)
elev_val       <- values(fire_elev.st)
elev2_val      <- values(fire_elev.st2)
lat_val        <- values(fire_lat.st)
size_val       <- values(fire_size)
cc_val         <- values(fire_cc.st)
precc_val      <- values(fire_precc.st)
precc_val[is.na(precc_val)] <- 0  # Code breaks with NA precc values. Replace with mean.
snag_cc_val    <- values(snag_fire_cc.st)
snag_cc2_val   <- values(snag_fire_cc.st2)
snag_precc_val <- values(snag_fire_precc.st)
snag_precc_val[is.na(snag_precc_val)] <- 0  # Code breaks with NA precc values. Replace with mean.
snag_size_val <- values(fire_size_raw)

# Start Pixel Loop
nsim <- 5000  # An arbitratily large number over which to sample from posterior
dens_val <- array(dim = c(length(pix), 2))  # Define output file

start.time <- Sys.time()
for(i in 1:length(pix)) {
  # Pick posterior values for this sim
  psi_i          <-round(runif(nsim, 1, length(post_psi[[2]])))
  distrib_psi    <-list(post_psi$b1[psi_i], 
                        post_psi$b2[psi_i], 
                        post_psi$b3[psi_i], 
                        post_psi$b4[psi_i], 
                        post_psi$b5[psi_i], 
                        post_psi$b6[psi_i])
  distrib_psi_b0 <- post_psi$b0[psi_i, ]  # intercepts based on 11 WHR-types
  hr_i           <- round(runif(nsim, 1, length(post_hr[[1]])))
  distrib_hr     <- list(post_hr$b0[hr_i], post_hr$b1[hr_i])
  snag_i         <- round(runif(nsim, 1, length(post_snag[[1]])))
  distrib_snag   <- list(post_snag$b0[snag_i],
                         post_snag$b1[snag_i],
                         post_snag$b2[snag_i],
                         post_snag$b3[snag_i],
                         post_snag$b4[snag_i],
                         post_snag$b5[snag_i],
                         post_snag$b6[snag_i])
  
  # Run model with chosen posterior
  snag_sim <- exp(distrib_snag[[1]] + 
                    distrib_snag[[2]] * snag_precc_val[pix[i]] + 
                    distrib_snag[[3]] * snag_cc_val[pix[i]] + 
                    distrib_snag[[4]] * snag_cc2_val[pix[i]] + 
                    distrib_snag[[5]] * snag_precc_val[pix[i]] * snag_cc_val[pix[i]] + 
                    distrib_snag[[6]] * snag_precc_val[pix[i]] * snag_cc2_val[pix[i]] + 
                    distrib_snag[[7]] * snag_size_val[pix[i]])
  snag_sim <- snag_sim * 0.2295687*10  # Convert to m2/ha
    
  psi_sim <- expit(distrib_psi_b0[, fire_veg2[pix[i]]] + 
                     distrib_psi[[1]] * elev_val[pix[i]] +  
                     distrib_psi[[2]] * elev2_val[pix[i]] + 
                     distrib_psi[[3]] * size_val[pix[i]] + 
                     distrib_psi[[4]] * cc_val[pix[i]] + 
                     distrib_psi[[5]] * precc_val[pix[i]] + 
                     distrib_psi[[6]] * lat_val[pix[i]])
  
  hr_sim               <- (exp(distrib_hr[[1]] + distrib_hr[[2]] * snag_sim))
  hr_sim[hr_sim < 20]  <- 20  # Set minimum home range as 20 ha
  hr_sim[hr_sim > 825] <- 825  # Set maximum home range as 825 ha
  hr_sim               <- 1 / hr_sim  # Density given occupancy, pairs per hectare
  hr_sim               <- hr_sim * (((xres(fire_cc)) ^ 2) / (100 ^ 2))  # Density given occupancy, pairs per pixel
  
  dens_vals <- hr_sim * psi_sim
  dens_val[i, ] <- c(exp(-1 * Mode(-log(dens_vals))), sd(dens_vals, na.rm = T))
  print(paste(round((i / length(pix)) * 100, 2), "%"))
}
end.time=Sys.time()
elapsed.time <- difftime(end.time, start.time, units='mins')
elapsed.time


#------------- Output values into rasters-------------
setwd(folder)
setwd("./Output")

val_frame <- rep(NA, length(in_fire))
raster_frame <- fire_outline_rast
values(raster_frame)[pix] <- 0

dens_mode_val      <- val_frame
dens_mode_val[pix] <- dens_val[, 1]
dens_mode          <- raster_frame
values(dens_mode)  <- dens_mode_val
dens_sd_val        <- val_frame
dens_sd_val[pix]   <- dens_val[, 2]
dens_sd            <- raster_frame
values(dens_sd)    <- dens_sd_val

writeRaster(dens_mode, "dens_mode.tif", format="GTiff", overwrite=T)
writeRaster(dens_sd, "dens_sd.tif", format="GTiff", overwrite=T)
