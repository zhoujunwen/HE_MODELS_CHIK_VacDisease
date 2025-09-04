####################################################################################################
# Example of code for analysis of chikungunya suitability, foi, incidence pipeline                 #
# Accompanying paper/preprint:  Health and economic burden of chikungunya infection and potential  #
# benefits of vaccination in 32 countries: a vaccine impact modelling study                        #
####################################################################################################

# for suitability code
library(devtools)
#install_github("SEEG-Oxford/seegSDM",lib="~/.local/rlibs")

library(seegSDM, lib.loc = "~/.local/rlibs")
library(raster)
library(rgdal)
library(sp)
library(gtools)
library(pivottabler)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rasterVis)
library(minpack.lm)

# for FOI code #
library("tidyverse")
library("maptools")
library("exactextractr")
library("terra")
library("sf")
library("geodata")
library(investr)
library(caret)
library(car)
library(parallel)
library(tibble)
library(purrr)
library(furrr)
library(future)
library(progressr)


setwd("~/Chik/CHIKV Modelling")

#load raster data
bb<-extent(-180,180,-60,85)
importrasters <- defmacro(layername,infile,
                          expr = {
                            layername <- raster(infile)
                            projection(layername) <- wgs84(TRUE)
                            layername<-setExtent(layername,bb,keepres=TRUE)
                          })



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

#import admin rasters, create a brick 
importrasters(admin0,"Admin/gau2015A0codes5km.tif")
importrasters(admin1,"Admin/gau2015A1codes5km.tif")
importrasters(admin2,"Admin/gau2015A2codes5km.tif")

admin<-brick(admin0,admin1,admin2)

#import admin area rasters
importrasters(area0,"Admin/areagadm0.tif")
importrasters(area1,"Admin/areagadm1.tif")
importrasters(area2,"Admin/areagadm2.tif")

area<-brick(area0,area1,area2)

#load occurrence data
library(readxl)

occurrence <- read_excel("CHIKV_Occurrence_Raw_2024Update.xlsx")

####TEMPORAL STANDARDISATION #####

#Extract GAUL Codes#

occurrence$Latitude<-as.numeric(occurrence$Latitude)
occurrence$Longitude<-as.numeric(occurrence$Longitude)
coordinates(occurrence)<- ~ Longitude + Latitude

GAULvals<-extract(admin,occurrence)
AREAvals<-extract(area,occurrence)

occ_withGAUL<-cbind(occurrence,GAULvals,AREAvals)

occ_withGAUL<-as.data.frame(occ_withGAUL)

occ_withGAUL$GAUL = ifelse(occ_withGAUL$Admin==0,occ_withGAUL$gau2015A0codes5km,
                           ifelse(occ_withGAUL$Admin==1,occ_withGAUL$gau2015A1codes5km,
                                  ifelse(occ_withGAUL$Admin==2,occ_withGAUL$gau2015A2codes5km,-999)))

occ_withGAUL$Area = ifelse(occ_withGAUL$Admin==0,occ_withGAUL$areagadm0,
                           ifelse(occ_withGAUL$Admin==1,occ_withGAUL$areagadm1,
                                  ifelse(occ_withGAUL$Admin==2,occ_withGAUL$areagadm2,-999)))

occ_withGAUL_clean<-subset(occ_withGAUL,occ_withGAUL$Area<2500)

#spdf <- SpatialPointsDataFrame(coords = xy, data = occ_withGAUL,
               #                proj4string = CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

# get date ranges (1 if same year)
range <- occ_withGAUL_clean$End_Year - occ_withGAUL_clean$Start_Year + 1

# make index for repetitions
n <- nrow(occ_withGAUL_clean)
rep_idx <- rep(1:n, times = range)

years <- occ_withGAUL_clean$Start_Year[rep_idx] + unlist(lapply(range, seq_len)) - 1

df <- occ_withGAUL_clean[rep_idx, ]
df <- df[, !(names(df) %in% c('Start_Year', 'End_Year'))]

df$Year <- as.numeric(years)

points<-subset(occ_withGAUL_clean,Admin==-999)
cellnums <- cellFromXY(admin, points[c('Longitude', 'Latitude')])
points <- cbind(cellnums, points)
points$duplicates <- duplicated(cbind(points$cellnums, points$Year))
points_std<-subset(points,duplicates==FALSE)
points_std <- subset(points_std, select=-c(cellnums,UNIQUEID,GAUL,Area,duplicates))
write.csv(points_std, file="Occurrence/points_std.csv")


polygons<-subset(occ_withGAUL_clean,Admin!=-999)
polygons$duplicates<-duplicated(cbind(polygons$GAUL,polygons$Year))
polygons_std<-subset(polygons,duplicates==FALSE)
points_std <- subset(polygons_std, select=-c(cellnums,UNIQUEID,GAUL,Area,duplicates))
write.csv(polygons_std, file="Occurrence/polygons_std.csv")

#convert occurrence data to spatial objects
polySPDF<-occurrence2SPDF(polygons_std)
pointsSPDF<-occurrence2SPDF(points_std)


#generate pseudo-absences and extract data

importrasters(aeg_unsuit,"bradtaegyptiunsuiteq1.tif")

set.seed(1)
#sample background points 
bg <- bgSample(aeg_unsuit,
              n = 10000,
              prob = TRUE,
              replace = TRUE,
              spatial = FALSE)
write.csv(bg, file="Occurrence/bg.csv")


colnames(bg) <- c('Longitude', 'Latitude')
bg <- data.frame(bg)


# get the covariate values
dat_points<-extract(covariates,pointsSPDF)
dat_pointS1 <- cbind(PA = rep(1, nrow(pointsSPDF)),
                     pointsSPDF@coords,
                     dat_points)

dat_poly<-extractAdmin(polySPDF,covariates,admin,fun='mean')
dat_poly2 <-cbind(PA = rep(1, nrow(polySPDF)),
                  polySPDF@coords,
                  dat_poly)

dat_bg <- extract(covariates, bg)
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
dat_all <- na.omit(dat)


##################################################################################################################
# run the ensemble 

nboots <- 100

# bootstrap subsample the dataset nboots times
data_list <- replicate(nboots,   
                       subsample(dat_all,  
                                 nrow(dat_all),  
                                 minimum = c(100, 100),  # make sure there are 100 presence/ 100 absence in each subsample
                                 prescol = 1,  
                                 replace = TRUE), 
                       simplify = FALSE) 

#set up a cluster 
sfInit(cpus=48,parallel=TRUE)
sfLibrary(seegSDM)

model_list<-sfLapply(data_list,
                     runBRT,
                     5:ncol(dat_all),
                     1,
                     covariates,
                     wt=2,
                     gbm.coords = 3:4,
                     var.monotone=c(1,1,1,1,1,1))

stat_lis <- sfLapply(model_list, getStats)

sfStop()


##################################################################################################################
# get stats 

stats <- do.call("rbind", stat_lis)
write.csv(stats, file="results/stats.csv")

boxplot(stats[, 3:7], col = "grey", ylim = c(0, 1))

relinf <- getRelInf(model_list, plot = TRUE)
par(mfrow = c(1, 3))
effect <- getEffectPlots(model_list, plot = TRUE)
write.csv(relinf,file="results/relinf.csv")


preds <- lapply(model_list, function(x) x$pred)
preds <- brick(preds)
preds_combined <- combinePreds(preds)

plot(preds_combined, zlim = c(0, 1))

preds_combined$uncertainty <- 
  preds_combined$quantile_0.975 - preds_combined$quantile_0.025


tiff("results/mean_suitability.tiff")
plot(preds_combined$mean, zlim = c(0, 1), main = "mean")
dev.off()

tiff("results/uncertainty_suitability.tiff")
plot(preds_combined$uncertainty, col = topo.colors(100), main = "uncertainty")
dev.off()

##################################################################################
###OUTPUT RESULTS###
writeRaster(preds_combined$mean, "results/mean_prediction", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.025, "results/lowerCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$quantile_0.975, "results/upperCI", format = "GTiff", overwrite=TRUE)
writeRaster(preds_combined$uncertainty, "results/uncertainty", format = "GTiff", overwrite=TRUE)

#####################################################################################################
# masking: 

bb<-extent(-180,180,-60,85)
importrasters <- defmacro(layername,
                          infile,
                          expr = {
                            layername <- raster(infile)
                            layername<-setExtent(layername,bb,keepres=TRUE)
                          })

importrasters(CHIKVsuit_unmasked,"~/Chik/CHIKV Modelling/results/mean_prediction.tif")
importrasters(tempsuit_aegypti,"BRADYTSI_AESUIT_AEGYPTI_2_5min_Rcp45_2015_clean.tif")
importrasters(tempsuit_albo,"BRADYTSI_AESUIT_ALBOPICTUS_2_5min_Rcp45_2015_clean.tif")
occurrence_std<-read.csv("occ_combined_std.csv")
coordinates(occurrence_std)= ~ x+y

## create a combined aedes temperature suitability map
aedes_combined_mask<- max(tempsuit_aegypti,tempsuit_albo) # have changed this to max() instead of just summing them
occ_tempsuit_vals<-extract(aedes_combined_mask,occurrence_std)

write.csv(occ_tempsuit_vals,file="occ_tempsuit_vals.csv")

## explore potential suitability thresholds based on % of occurrences that fall within
quantile(occ_tempsuit_vals, probs = c(.01,.025,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50), na.rm=TRUE)
threshold_90 <- quantile(occ_tempsuit_vals, probs = 0.10, na.rm=TRUE)
# threshold is 0.4886158 (so approximately 0.5)

## this masks according to a set threshold - change the value to coincide with which threshold is chosen 
tempsuit_mask_threshold<-aedes_combined_mask>threshold_90
CHIKVsuit_mask_threshold<-tempsuit_mask_threshold*CHIKVsuit_unmasked

writeRaster(CHIKVsuit_mask_threshold, "results/CHIKVsuit_2500_mask_threshold90.tif", format = "GTiff", overwrite=TRUE)


####################################################################################################
# save preds to be able to call without rerunning if it does crash:

saveRDS(preds, "results/suitability_preds.RDS")
saveRDS(shp_country, "results/shp_country.RDS")

rm(list=setdiff(ls(), c("preds", "shp_country","threshold_90")))

preds <- readRDS("results/suitability_preds.RDS")
shp_country <- readRDS("results/shp_country.RDS")


######################################### FOI estimates ########################################################

wd <- getwd()
path_data <- file.path("/data/Chik")
path_output <- file.path("~/Chik/output")
path_tmp <- file.path(path_data, "tmp")

# set raster and future options
rasterOptions(tmpdir = path_tmp)
options(future.globals.maxSize = 25 * 1024^3) 

extract <- raster::extract
select <- dplyr::select

# relevant countries:
nsuit <- 100
nsim <- 1
ndraw <- 100
ncores <- 15



# number of countries: 
#shp_country <- getData('countries')
ncountries <- length(unique(shp_country@data$ISO))

incidence_result <- data.frame(
  country = rep(levels(shp_country@data$ISO), nsim * ndraw * nsuit),
  nsuit = rep(1:nsuit, each = nsim * ndraw * ncountries),
  draw = rep(rep(rep(1:ndraw, each = ncountries), each = nsim), nsuit),
  sim = rep(rep(rep(1:nsim, each = ncountries), ndraw), nsuit),
  incidence = rep(rep(NA, ncountries * nsim * ndraw), nsuit)
)

# NLS
opt_dat <- "FULL" # specify if excluding high suitability very low foi as PART
# otherwise PART


foi_all <- readRDS("/data/Chik/foi_all_2000.RDS")


rst_hpdr <- raster(file.path(path_data, "ppp_2019_1km_Aggregated.tif"))
# change to 5km:
rst_hpdr <- terra::aggregate(rst_hpdr, 5, fun = sum)


control_list <- nls.control(maxiter = 50000)


# Process combination function with detailed debug messages
process_combination <- function(h, i, shp_country, rst_hpdr, control_list, foi_all, path_data,
                                 preds,threshold_90) {
  message(sprintf("Processing suit = %d, draw = %d", h, i))
  options(future.globals.maxSize = 25 * 1024^3)
  foi_all_1 <- foi_all %>% filter(draw %in% c(i))
  
  # given working in a parrallel process, make unique saves for each one
  temp_file <- tempfile(paste0("pred_", i), fileext = ".tif")
  on.exit(unlink(temp_file), add = TRUE)  # Ensures file deletion after function execution
  writeRaster(preds[[i]], filename = temp_file, format="GTiff", overwrite = TRUE)
  sd_suit <- raster(temp_file)
  crs(sd_suit) <- CRS("+proj=longlat +datum=WGS84 +no_defs") 
  
  sd_suit[sd_suit < threshold_90] <- 0
  
  foi_all_1 <- st_as_sf(foi_all_1)
  foi_all_1 <- as(foi_all_1, "Spatial")
  # # reprojecting preds takes ages, so probably do something with foi_all_1 instead:
  foi_all_1 <- spTransform(foi_all_1, CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  tmp1 <- exact_extract(sd_suit, foi_all_1)
  if (length(tmp1) == 0) stop("Exact extract returned no data for suitability")
  tmp2 <- map(tmp1, function(x) {
    temp1 <- x$value[!is.na(x$value)]
    if (length(temp1) > 0) {
      output <- tibble(val = mean(temp1, na.rm = TRUE))
    } else {
      output <- tibble(val = NA)
    }
    return(output)
  })
  message("Step 1 complete")
  
  temp1 <- exact_extract(rst_hpdr, foi_all_1)
  if (length(temp1) == 0) stop("Exact extract returned no data for population")
  temp2 <- map(temp1, function(x) {
    temp <- x$value[!is.na(x$value)]
    if (length(temp) > 0) {
      output <- tibble(pop = sum(temp, na.rm = TRUE))
    } else {
      output <- tibble(pop = NA)
    }
    return(output)
  })
  message("Step 2 complete")
  

  # Calculate population-weighted suitability by country

  # Ensure country polygons match CRS
  shp_country <- spTransform(shp_country, CRS(projection(sd_suit)))
  
  # Resample population raster to match suitability raster 
  rst_hpdr_resampled <- resample(rst_hpdr, sd_suit, method = "bilinear")
  
  # Multiply: suitability * population
  suitability_weighted <- overlay(sd_suit, rst_hpdr_resampled, fun = function(suit, pop) suit * pop)
  
  rstack_country <- stack(suitability_weighted, rst_hpdr_resampled)
  
  # Extract population-weighted suitability per country polygon
  country_suitability <- exact_extract(rstack_country, shp_country, 
                                       fun = function(df) {
                                         total_pop <- sum(df[[2]], na.rm = TRUE)
                                         if (total_pop == 0) return(NA)
                                         weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                                         return(weighted_mean)
                                       },
                                       summarize_df = TRUE)
  
  country_results <- tibble(
    country = shp_country@data$NAME_ENGLISH,
    pop_weighted_suitability = country_suitability,
    suit_draw = i  # optional: track draw or bootstrap iteration
  )

  tmp3 <- map(1:length(tmp2),
              ~bind_cols(as_tibble(foi_all_1@data)[.x,], tmp2[[.x]], temp2[[.x]])
  ) %>%
    bind_rows() %>%
    filter(!is.na(val)) %>%
    mutate(pop = replace_na(pop, 0)) %>%
    group_by(study_id, COUNTRY, NAME_1, NAME_2, NAME_3, draw, foi2000) %>%
    mutate(p = pop / sum(pop)) %>%
    summarize(val = sum(val * p), pop = sum(pop), .groups = "drop")
  message("Step 3 complete")
  
  tp2 <- tmp3 %>%
    rename(suit = val) %>%
    group_by(study_id, COUNTRY, NAME_1, NAME_2, NAME_3, draw, foi2000) %>%
    mutate(p = pop / sum(pop)) %>%
    summarize(suit = sum(suit * p),
              pop = sum(pop), .groups = "drop")
  message("Step 4 complete")
  
  tp3 <- tp2 %>%
    filter(!is.na(suit)) %>%
    group_by(study_id, COUNTRY, NAME_1, NAME_2, NAME_3, draw) %>%
    mutate(p = pop / sum(pop)) %>%
    summarise(foi = sum(foi2000 * p),
              suit = sum(suit * p))
  message("Step 5 complete")
  
  tmp <- tp3 %>%
    ungroup %>%
    filter(foi == 0) %>%
    sample_frac(size = 0.2, replace = FALSE) %>%
    dplyr::select(study_id, COUNTRY, NAME_1, NAME_2, NAME_3, draw) %>%
    mutate(ex = 1)
  df <- tp3 %>%
    left_join(tmp, by = c("study_id", "COUNTRY", "NAME_1", "NAME_2", "NAME_3", "draw")) %>%
    filter(is.na(ex)) %>%
    rename(y = foi, x = suit) %>%
    ungroup()
  message("Step 6 complete")
  
  df2 <- df #%>% filter(!(x > 0.75 & y < 0.005))
  df2$y <- ifelse(df2$x <as.numeric(threshold_90), 0, df2$y)
  # truncate at 97.5%
  truncate975 <- as.numeric(quantile(df2$y, probs = 0.975))
  df2$y <- ifelse(df2$y > truncate975, truncate975, df2$y)
  
  
  start_values <- list(Asym = max(df2$y)/2,
                       xmid = median(df2$x[df2$y!=0]),
                       scal = diff(range(df2$x))/3)
  
  model <- tryCatch({
    nlsLM(y ~ SSlogis(x, Asym, xmid, scal),
          data = df2,
          start = start_values,
          lower = c(Asym=0, xmid = min(df2$x), scal=0.01),
          upper = c(Asym=0.04, xmid = max(df2$x), scal=0.5))
    
  }, error = function(e) NULL)

  
  if (is.null(model)) {
    message(sprintf("Model fitting failed for h = %d, i = %d", h, i))
    return(NULL)
  }
  message("Model fitting complete")
  
  suit_vals <- values(sd_suit)
  # Mask out suitability below threshold or NA
  valid_mask <- !is.na(suit_vals) & suit_vals >= threshold_90
  # Make prediction input:
  b <- tibble(x = suit_vals)
  c <- rep(0, length(suit_vals))  # initialize as 0
  c[valid_mask] <- predict(model, newdata = b[valid_mask, , drop = FALSE])
  pred <- as.vector(c)
  
  pred[suit_vals == 0 | is.na(suit_vals)] <- 0 # this works because a and c in the same order

  message("Prediction complete")
  
  #######################################
  process_simulation_extended <- function(h, i, pred, sd_suit, rst_hpdr, shp_country) {
    try({
      message(sprintf("Processing draw %d and suit %d", h, i))
      
      rst_foiPred <- sd_suit
      values(rst_foiPred) <- pred
      
      rst_foiPred_hpdr <- resample(rst_foiPred, rst_hpdr)
      
      rst_incidence_hpdr <- overlay(rst_foiPred_hpdr, rst_hpdr, fun = function(x, y) {
        return(x * y)
      })
      
      dat_incidence_by_country <- exact_extract(rst_incidence_hpdr, shp_country)
      
      df_incidence_by_country <- do.call(rbind, lapply(seq_along(dat_incidence_by_country), function(idx) {
        x <- dat_incidence_by_country[[idx]]
        country_iso <- shp_country@data$ISO[idx]
        
        valid_values <- x$value[!is.na(x$value)]
        
        if (length(valid_values) > 0) {
          output <- tibble(
            id = 1:length(valid_values),
            incidence = valid_values,
            country = country_iso
          )
        } else {
          output <- tibble(
            id = 1,
            incidence = 0,
            country = country_iso
          )
        }
        
        return(output)
      }))
      
      summarized_data <- df_incidence_by_country %>%
        mutate(draw = h, suit = i) %>%  
        group_by(country, draw, suit) %>%
        summarize(incidence = sum(incidence, na.rm = TRUE), .groups = 'drop')
      
      return(summarized_data)
    }, silent = FALSE)
  }
  
  all_extended_results <- list()
  
  result <- process_simulation_extended(h, i, pred, sd_suit, rst_hpdr, shp_country)
  all_extended_results[[length(all_extended_results) + 1]] <- result
  
  final_extended_incidence_data <- bind_rows(all_extended_results)
  
  
  out <- list(h = h, i = i, results = final_extended_incidence_data)
  return(out)
}


##############################################################

cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(tidyverse)
  library(raster)
  library(exactextractr)
  library(terra)
  library(sf)
  library(investr)
  library(caret)
  library(car)
  library(tibble)
  library(dplyr)
  library(purrr)
  library(furrr)
  library(future)
  library(progressr)
  library(minpack.lm)
  options(future.globals.maxSize = 25 * 1024^3)
})


clusterExport(cl, c("foi_all", "shp_country", "rst_hpdr", "nsuit", "ndraw", "path_data", "control_list", "process_combination",
                     "preds","threshold_90"))

# Running the process in parallel
results <- parLapply(cl, 1:nsuit, function(h) {
  lapply(1:ndraw, function(i) {
    process_combination(h, i,  shp_country, rst_hpdr, control_list, foi_all, path_data, 
                         preds,threshold_90)
  })
})

stopCluster(cl)

saveRDS(results, file = "results/incidence_country_full_masked90_foi.RDS")
results <- readRDS("results/incidence_country_full_masked90_foi.RDS")

str(results)

check <- results[[2]][[1]][["results"]]
check%>%filter(country=="IND")
