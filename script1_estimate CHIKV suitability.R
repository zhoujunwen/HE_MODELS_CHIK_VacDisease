####################################################################################################
# Example of code for analysis of chikungunya suitability, foi, incidence pipeline                 #
# Accompanying paper/preprint:  Health-economic burden of chikungunya infection and potential      #
# cost-effectiveness of preventive vaccination in 31 countries                                     #
####################################################################################################

# START ----

# Clear Environment
rm(list=ls())

# Set work directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # the directory of this script

# load required packages

# Install required packages
# install.packages("terra")
# install.packages("devtools")
# library(devtools)
# Need to save the seegSDM packages somewhere first
# remotes::install_local("C:\\Users\\junwenz\\Documents\\seegSDM_0.1-9\\seegSDM")

# for suitability code
library(tidyverse)
library(terra)
library(sf)
library(seegSDM)
library(gtools)
library(parallel)
library(furrr)
library(future)
library(progressr)
# library(car)
# library(caret)

# Load data ----

# Raster - Covariate
bb<-extent(-180,180,-60,85)
importrasters <- defmacro(layername,infile,
                          expr = {
                            layername <- raster(infile)
                            projection(layername) <- wgs84(TRUE)
                            layername<-setExtent(layername,bb,keepres=TRUE)
                          })
setwd(wd_data)
#import covariate rasters, create a brick and plot the rasters
importrasters(evi,"Covariates/wg1115mx.tif")
importrasters(lst,"Covariates/wg1107mn.tif")
#band_1.1
importrasters(precip,"Covariates/wd1920a05.tif")
importrasters(pop,"Covariates/wdwpop5km.tif")
#band_1.2
importrasters(aedes,"Covariates/E4AEGUNS4BRMEAN15KbalV2julallvarmeanrfbrt.tif")
#band_1.3
importrasters(albo,"Covariates/ed410xyalboedafbal10meanbrtrfjul23.tif")
covariates<-brick(evi,lst,precip,pop,aedes,albo)

# Raster - Admin
importrasters(admin0,"Admin/gau2015A0codes5km.tif")
importrasters(admin1,"Admin/gau2015A1codes5km.tif")
importrasters(admin2,"Admin/gau2015A2codes5km.tif")
admin<-brick(admin0,admin1,admin2)

# Raster - Area
importrasters(area0,"Admin/areagadm0.tif")
importrasters(area1,"Admin/areagadm1.tif")
importrasters(area2,"Admin/areagadm2.tif")
area<-brick(area0,area1,area2)

# Occurrence data
library(readxl)
occurrence <- read_excel("CHIKV_Occurrence_Raw_2024Update.xlsx")

# Temperature suitability for aegypti and albo
# read in temperature suitability and get rid of those that have a temp suitability of 0: 
## create a combined aedes temperature suitability map
importrasters(tempsuit_aegypti,"Tempsuit/BRADYTSI_AESUIT_AEGYPTI_2_5min_Rcp45_2015_clean.tif")
importrasters(tempsuit_albo,"Tempsuit/BRADYTSI_AESUIT_ALBOPICTUS_2_5min_Rcp45_2015_clean.tif")
aedes_combined_mask<- max(tempsuit_aegypti,tempsuit_albo) # have changed this to max() instead of just summing them

# Clean occurrence data ----

# Extract GAUL Codes#
occurrence$Latitude<-as.numeric(occurrence$Latitude)
occurrence$Longitude<-as.numeric(occurrence$Longitude)
coordinates(occurrence)<- ~ Longitude + Latitude

# extract combined temperature suitability at each occurrence point
ts_vals <- terra::extract(aedes_combined_mask, occurrence)
occurrence$tempsuit <- ts_vals
GAULvals<-terra::extract(admin,occurrence)
AREAvals<-terra::extract(area,occurrence)
occ_withGAUL<-cbind(occurrence,GAULvals,AREAvals)
occ_withGAUL<-as.data.frame(occ_withGAUL)
occ_withGAUL$GAUL = ifelse(occ_withGAUL$Admin==0,occ_withGAUL$gau2015A0codes5km,
                           ifelse(occ_withGAUL$Admin==1,occ_withGAUL$gau2015A1codes5km,
                                  ifelse(occ_withGAUL$Admin==2,occ_withGAUL$gau2015A2codes5km,-999)))
occ_withGAUL$Area = ifelse(occ_withGAUL$Admin==0,occ_withGAUL$areagadm0,
                           ifelse(occ_withGAUL$Admin==1,occ_withGAUL$areagadm1,
                                  ifelse(occ_withGAUL$Admin==2,occ_withGAUL$areagadm2,-999)))
occ_withGAUL_clean<-subset(occ_withGAUL,occ_withGAUL$Area<2500 & tempsuit > 0)
# add in here that temperature suitability should be above 0. 

# get date ranges (1 if same year)
range <- occ_withGAUL_clean$End_Year - occ_withGAUL_clean$Start_Year + 1

# make index for repetitions (rep each index 'range' times)
n <- nrow(occ_withGAUL_clean)
rep_idx <- rep(1:n, times = range)

# calculate years (start_year + 0:range)
years <- occ_withGAUL_clean$Start_Year[rep_idx] + unlist(lapply(range, seq_len)) - 1

# subset dataframe
df <- occ_withGAUL_clean[rep_idx, ]

# remove start and end year columns
df <- df[, !(names(df) %in% c('Start_Year', 'End_Year'))]

# add years
df$Year <- as.numeric(years)

# Clean point
points<-subset(occ_withGAUL_clean,Admin==-999)
cellnums <- cellFromXY(admin, points[c('Longitude', 'Latitude')])
points <- cbind(cellnums, points)
points$duplicates <- duplicated(cbind(points$cellnums, points$Year))
points_std<-subset(points,duplicates==FALSE)
points_std <- subset(points_std, select=-c(cellnums,UNIQUEID,GAUL,Area,duplicates))


# Clean polygons
polygons<-subset(occ_withGAUL_clean,Admin!=-999)
polygons$duplicates<-duplicated(cbind(polygons$GAUL,polygons$Year))
polygons_std<-subset(polygons,duplicates==FALSE)
# polygons_std <- subset(polygons_std, select=-c(UNIQUEID,GAUL,Area,duplicates))

#convert occurrence data to spatial objects
polySPDF<-occurrence2SPDF(polygons_std)
pointsSPDF<-occurrence2SPDF(points_std)

setwd(wd_output)
write.csv(points_std, file="points_std.csv")
write.csv(polygons_std, file="polygons_std.csv")

# Generate pseudo-absences ----

#generate pseudo-absences and extract data
setwd(wd_data)
importrasters(aeg_unsuit,"Tempsuit/bradtaegyptiunsuiteq1.tif")

set.seed(1)
#sample background points 
bg <- bgSample(aeg_unsuit,
              n = 10000,
              prob = TRUE,
              replace = TRUE,
              spatial = FALSE)
setwd(wd_output)
write.csv(bg, file="bg.csv")
# bg <- read.csv(file = "bg.csv")[,-1]
colnames(bg) <- c('Longitude', 'Latitude')
bg <- data.frame(bg)

# get the covariate values
dat_points<-terra::extract(covariates,pointsSPDF)
dat_pointS1 <- cbind(PA = rep(1, nrow(pointsSPDF)),
                     pointsSPDF@coords,
                     dat_points)

dat_poly<-extractAdmin(polySPDF,covariates,admin,fun='mean')
dat_poly2 <-cbind(PA = rep(1, nrow(polySPDF)),
                  polySPDF@coords,
                  dat_poly)

dat_bg <- terra::extract(covariates, bg)
dat_bg2 <-cbind(PA = rep(0, nrow(bg)),
                bg,
                dat_bg)

occ<-rbind(dat_pointS1,dat_poly2)

#tell it that dat and BG have equal weights
occ <-cbind(occ,weights= rep(1, nrow(occ)))
bg2 <- as.data.frame(dat_bg2)
weights <- (nrow(occ)/nrow(bg2))
bg3<-cbind(bg2,weights)

#combine the occurrence and background records
dat <-rbind(occ,bg3)
dat<- dat[,c(1,10,2,3,4,5,6,7,8,9)]
dat<-as.data.frame(dat)

# remove NAs
dat_all <- na.omit(dat)

# RUN MODELS ----

# define number of bootstraps
nboots <- 100

# bootstrap subsample the dataset nboots times
data_list <- replicate(nboots,   # replicate this procedure nboots times
                       subsample(dat_all,  # the full dataset
                                 nrow(dat_all),  # the number of rows in the dataset
                                 minimum = c(100, 100),  # make sure there are 100 presence/ 100 absence in each subsample
                                 prescol = 1,  # which column number gives the presence/absence column
                                 replace = TRUE),  # need to sample with replacement for true bootstrap
                       simplify = FALSE)  # simplify as a list

#set up a cluster of cpus in with parallel execution
sfInit(cpus=30,parallel=TRUE)
sfLibrary(seegSDM)

#fit the models of the ensemble in parallel
model_list<-sfLapply(data_list,
                     runBRT,
                     5:ncol(dat_all),
                     1,
                     covariates,
                     wt=2,
                     gbm.coords = 3:4,
                     var.monotone=c(1,1,1,1,1,1))

# get cv statistics in parallel
stat_lis <- sfLapply(model_list, getStats)

# Now we've finished with the parallel cluster, we should shut it down
sfStop()

# Output ----

# SUMMARIZE THE ENSEMBLE

# convert the list into a matrix using the do.call function
stats <- do.call("rbind", stat_lis)

setwd(wd_output)
write.csv(stats, file="stats.csv")
# look at it
head(stats)

# and produce a boxplot of a few imnportant statistics
boxplot(stats[, 3:7], col = "grey", ylim = c(0, 1))

relinf <- getRelInf(model_list, plot = TRUE)
par(mfrow = c(1, 3))
effect <- getEffectPlots(model_list, plot = TRUE)
write.csv(relinf,file="relinf.csv")

# lapply to extract the predictions into a list
preds <- lapply(model_list, function(x) x$pred)
saveRDS(preds, "suitability_preds_mask0.RDS")

preds <- readRDS(file.path(wd_output, "suitability_preds_mask0.RDS"))

# coerce the list into a rasterbrick
preds <- brick(preds)

# now we can run combinePreds
preds_combined <- combinePreds(preds)

# plot the resulting maps
plot(preds_combined, zlim = c(0, 1))

# calculate uncertainty
preds_combined$uncertainty <- #preds[[4]] - preds[[3]]
  preds_combined$quantile_0.975 - preds_combined$quantile_0.025

# plot mean and uncertainty
tiff("mean_suitability.tiff")
# plot mean
plot(preds_combined$mean, zlim = c(0, 1), main = "mean")
dev.off()
# and uncertainty
tiff("uncertainty_suitability.tiff")
plot(preds_combined$uncertainty, col = topo.colors(100), main = "uncertainty")
dev.off()

# OUTPUT RESULTS 
writeRaster(preds_combined$mean, "mask0_mean", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.025, "mask0_lowerCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.975, "mask0_upperCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$uncertainty, "mask0_uncertainty", format = "GTiff", overwrite=TRUE)

# MASKING ----

# masking: 
preds <- readRDS("suitability_preds_mask0.RDS")
#load raster data and standardise extent
bb<-extent(-180,180,-60,85)
importrasters <- defmacro(layername,
                          infile,
                          expr = {
                            layername <- raster(infile)
                            layername<-setExtent(layername,bb,keepres=TRUE)
                          })
# importrasters(CHIKVsuit_unmasked,"results/mean_prediction.tif")
importrasters(tempsuit_aegypti,file.path(wd_data, "Tempsuit/BRADYTSI_AESUIT_AEGYPTI_2_5min_Rcp45_2015_clean.tif"))
importrasters(tempsuit_albo,file.path(wd_data, "Tempsuit/BRADYTSI_AESUIT_ALBOPICTUS_2_5min_Rcp45_2015_clean.tif"))
aedes_combined_mask<- max(tempsuit_aegypti,tempsuit_albo) # have changed this to max() instead of just summing them

aedes_fl_nomask <- aedes_combined_mask > 0

# For each prediction, mask temp suit = 0
preds2 <- map(preds, ~.x * aedes_fl_nomask)

saveRDS(preds2, "suitability_preds_mask1.RDS")

# coerce the list into a rasterbrick
preds2 <- brick(preds2)

# now we can run combinePreds
preds_combined <- combinePreds(preds2)

# plot the resulting maps
plot(preds_combined, zlim = c(0, 1))

# calculate uncertainty
preds_combined$uncertainty <- #preds[[4]] - preds[[3]]
  preds_combined$quantile_0.975 - preds_combined$quantile_0.025

# plot mean and uncertainty
tiff("mean_suitability.tiff")
# plot mean
plot(preds_combined$mean, zlim = c(0, 1), main = "mean")
dev.off()
# and uncertainty
tiff("uncertainty_suitability.tiff")
plot(preds_combined$uncertainty, col = topo.colors(100), main = "uncertainty")
dev.off()

# OUTPUT RESULTS 
writeRaster(preds_combined$mean, "mask1_mean", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.025, "mask1_lowerCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.975, "mask1_upperCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$uncertainty, "mask1_uncertainty", format = "GTiff", overwrite=TRUE)
