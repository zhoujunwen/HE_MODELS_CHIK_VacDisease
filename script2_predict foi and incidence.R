####################################################################################################
# Example of code for analysis of chikungunya suitability, foi, incidence pipeline                 #
# Accompanying paper/preprint:  Health-economic burden of chikungunya infection and potential      #
# cost-effectiveness of preventive vaccination in 31 countries                                     #
####################################################################################################

# Start ----

# Clear Environment
rm(list=ls())

# Set work directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # the directory of this script

# Package
library(tidyverse)
library(terra)
library(sf)
library(raster)
library(exactextractr) # Extraction from Raster
library(car)
library(minpack.lm)
# library("parallel")
# library("doRNG")
# library("doParallel")
select <- dplyr::select

# Prepare data ----

## FOI ----
  
foi_est_para <- readRDS(file.path(wd_data, "foi_all.rds")) %>% 
  select(loc_id = Location_uni_id, para1, para2)

foi_est <- list()
set.seed(1234)
for(i in 1:nrow(foi_est_para)){
  temp <- foi_est_para[i,]
  foi_est[[i]] <- tibble(loc_id = temp$loc_id,
                         sim = 1:100,
                         foi = rbeta(100, temp$para1, temp$para2))
}
foi_est <- bind_rows(foi_est) %>% arrange(sim, loc_id) %>% relocate(sim, loc_id)

## Location of FOI ----

# Foi location raster (data in tp)
loc_foi <- readRDS(file.path(wd_data, "foi_all.rds")) %>% 
  select(loc_id = Location_uni_id, geometry) %>% 
  arrange(loc_id)
loc_foi <- st_as_sf(loc_foi)
loc_foi <- as(loc_foi, "Spatial")
# # reprojecting preds takes ages, so probably do something with foi_all_1 instead:
loc_foi <- spTransform(loc_foi, CRS("+proj=longlat +datum=WGS84 +no_defs"))

## Location of FOI KG ----

loc_foi2 <- readRDS(file.path(wd_data, "foi_KG_loc.rds")) %>% 
  select(loc_id, geometry) %>% 
  arrange(loc_id)
loc_foi2 <- st_as_sf(loc_foi2)
loc_foi2 <- as(loc_foi2, "Spatial")
# # reprojecting preds takes ages, so probably do something with foi_all_1 instead:
loc_foi2 <- spTransform(loc_foi2, CRS("+proj=longlat +datum=WGS84 +no_defs"))

## Population density ----

# Population density in 5km 
hpdr <- raster(file.path(wd_data, "human_population_density_2019CRS.tif"))
hpdr <- terra::aggregate(hpdr, 5, fun = sum) # Have to do it like this and can not save the 5km hpdr then reload

## Country of interest ----

# shp_country <- readRDS(file.path(wd_data, "shp_country.RDS"))
# temp0 <- read.csv(file.path(wd_data, "final_countries.csv")) %>% filter(!GID_0 %in% c("TON", "MDV"))
# shp_country <- shp_country[shp_country@data$ISO %in% temp0$GID_0, ]
# saveRDS(shp_country, file.path(wd_data, "shp_country_final.RDS")) 
shp_country <- readRDS(file.path(wd_data, "shp_country_final.RDS"))

# Population density in the country of interest
pop_at_country <- exact_extract(hpdr, shp_country, 
                                fun = function(df) {return(sum(df, na.rm = TRUE))},
                                summarize_df = TRUE)

## Suitability data ----

# 100 predicted suitability data
suit_est <- readRDS(file.path(wd_data, "suitability_preds_mask1.RDS"))

# 1 sensitivty suitability data from Lim 2025 to be used below

# Mapping + prediction ----

tp <- list()
test <- list()
# For each suitability data
for(i in 1:101){
  
  if(i < 101){
    # Project suitability from unit = meter to unit = degree
    temp_file <- tempfile(paste0("suit_est", i), fileext = ".tif")
    on.exit(unlink(temp_file), add = TRUE)  # Ensures file deletion after function execution
    writeRaster(suit_est[[i]], filename = temp_file, format="GTiff", overwrite = TRUE)
    suit_i <- raster(temp_file)
    crs(suit_i) <- CRS("+proj=longlat +datum=WGS84 +no_defs") 
  } else {
    suit_i <- raster(file.path(wd_data, "CHIK_riskmap_wmean_masked.tif"))
    crs(suit_i) <- CRS("+proj=longlat +datum=WGS84 +no_defs") 
  }
  
  # Population weighted suitability at the country
  
  # > Make suitability the same resolutation as population density
  suit_i_resampled <- resample(suit_i, hpdr)
  
  # > Assign weight (population density) to suitability
  rstack_suit_i_weighted <- overlay(suit_i_resampled, hpdr, fun = function(suit, pop) {return(suit * pop)})
  
  # > Stack weighted suitability and total weight
  rstack_suit_hpdr <- stack(rstack_suit_i_weighted, hpdr)
  
  # > Derive the population weighted suitability by using weighted suitability dividing total weight
  suit_i_at_country <- exact_extract(rstack_suit_hpdr, shp_country, 
                                     fun = function(df) {
                                       total_pop <- sum(df[[2]], na.rm = TRUE)
                                       if (total_pop == 0) return(NA)
                                       weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                                       return(weighted_mean)
                                       },
                                       summarize_df = TRUE)
  temp0 <- exact_extract(suit_i, loc_foi2)
  temp <- map_dbl(temp0, ~mean(replace_na(.x$value, 0)))
  temp2 <- exact_extract(rstack_suit_hpdr, loc_foi2, 
                         fun = function(df) {
                           total_pop <- sum(df[[2]], na.rm = TRUE)
                           if (total_pop == 0) return(NA)
                           weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                           return(weighted_mean)
                           },
                         summarize_df = TRUE)
  
  test[[i]] <- tibble(loc_id = loc_foi2@data$loc_id,
                      suit_i_mean = temp,
                      suit_i_popweightd = temp2)
  
  # Mapping suitability to FOI
  
  # > Derive the suitability at FOI location
  temp <- exact_extract(suit_i, loc_foi)
  suit_i_at_loc <- map_df(temp1, ~tibble(suit = mean(replace_na(.x$value, 0))), .id = "loc_id") %>%
    mutate(loc_id = as.numeric(loc_id))
  temp <- exact_extract(rstack_suit_hpdr, loc_foi,
                        fun = function(df) {
                          total_pop <- sum(df[[2]], na.rm = TRUE)
                          if (total_pop == 0) return(NA)
                          weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                          return(weighted_mean)
                        },
                        summarize_df = TRUE)
  suit_i_at_loc2 <- tibble(loc_id = as.numeric(loc_foi@data$loc_id), suit = temp)

  # # > map suitability to each simulated FOI

  # df_output <- list()
  # for(j in 1:100){

  n_cores <- 50
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ptm <- proc.time()
  df_output <- foreach(
    j = 1:100,
    # i = 1:1000,
    # .combine = rbind,
    .packages = c("tidyverse", "terra", "sf", "raster", "exactextractr", "car", "minpack.lm")
    ) %dorng% {

      foi_j <- foi_est %>% filter(sim == j) %>% select(-sim)
      tp <- list()
      for(k in 1:2){

        temp <- list(suit_i_at_loc, suit_i_at_loc2)[[k]]
        suit_i_foi_j <- foi_j %>%
          left_join(temp, by = "loc_id") %>%
          select(loc_id, suit, foi)
        truncate975 <- as.numeric(quantile(suit_i_foi_j$suit, probs = 0.975))
        df <- suit_i_foi_j %>%
          mutate(y = foi,
                 x = ifelse(suit > truncate975, truncate975, suit)) %>%
          select(y, x)
        # start_values <- list(Asym = 0.03590443, xmid = 0.3217752, scal = 0.09638236)
        start_values <- list(Asym = max(df$y)/2,
                             xmid = median(df$x[df$y!=0]),
                             scal = diff(range(df$x))/3)
        control_list <- nls.control(maxiter = 50000)
        fl_bkmod <- 0
        model <- tryCatch({
          nls(y ~ SSlogis(x, Asym, xmid, scal), data = df, control = control_list, algorithm = "port",
              start = start_values,
              lower = c(Asym=0, xmid = min(df$x), scal =0.01), # set sensible lower bounds,
              # asymptote cannot be negative
              # xmid withint the range of x
              # positive scale (determines the steepness, which has to be positive)
              upper = c(Asym=max(df$y), xmid = max(df$x), scal=0.5)
              # asym more bounded so I don't get an annual force of infection (average over long period) above 0.05
              # xmid within the range of x
              # scale less then than range of x
          )
        }, error = function(e) NULL)

        if (is.null(model)) {
          message(sprintf("Model fitting failed for i = %d, j = %d", i, j))
          model <- tryCatch({
            nlsLM(y ~ SSlogis(x, Asym, xmid, scal),
                  data = df,
                  start = start_values,
                  lower = c(Asym=0, xmid = min(df$x), scal=0.01),
                  upper = c(Asym=0.04, xmid = max(df$x), scal=0.5))
          }, error = function(e) NULL)
          fl_bkmod <- 1
          return(NULL)
        }

        if (is.null(model)) {
          message(sprintf("Model 2 fitting failed for i = %d, j = %d", i, j))
          fl_bkmod <- 2
          return(NULL)
        }

        # Predict FOI based on the suitability map

        # > Extract all suitability values
        suit_vals <- values(suit_i)

        # > Mask out suitability NA (or above threshold if needed)
        # valid_mask <- !is.na(suit_vals) & suit_vals >= threshold_90
        valid_mask <- !is.na(suit_vals)

        # > Make prediction input:
        b <- tibble(x = suit_vals)
        c <- rep(0, length(suit_vals))  # initialize as 0

        # > Make prediction
        c[valid_mask] <- predict(model, newdata = b[valid_mask, , drop = FALSE])
        pred <- as.vector(c)
        pred[suit_vals == 0 | is.na(suit_vals)] <- 0 # this works because a and c in the same order

        # > Assign prediction to the map
        pred_i_j <- suit_i # FOI same map as suit map
        values(pred_i_j) <- pred # Assign predicted FOI to the corresponding suit point in the map

        # > Resample the FOI map to be the same resolution as population density
        pred_i_j_hpdr <- resample(pred_i_j, hpdr)

        # Extract country weighted FOI

        # > Assign weight to FOI (Population density)
        rstack_foi_i_j_weighted <- overlay(pred_i_j_hpdr, hpdr, fun = function(foi, pop){return(foi * pop)})

        # > Stack weighted FOI and total weight
        rstack_foi_i_j_hpdr <- stack(rstack_foi_i_j_weighted, hpdr)

        # > Derive the predicted population weighted FOI across country of interest using weighted suitability dividing total weight
        foi_i_j_at_country <- exact_extract(rstack_foi_i_j_hpdr, shp_country,
                                            fun = function(df) {
                                              total_pop <- sum(df[[2]], na.rm = TRUE)
                                              if (total_pop == 0) return(NA)
                                              weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                                              return(weighted_mean)
                                            },
                                            summarize_df = TRUE)

        # > Derive the predicted population weighted FOI across location of interest using weighted suitability dividing total weight
        foi_i_j_at_loc <- exact_extract(rstack_foi_i_j_hpdr, loc_foi,
                                            fun = function(df) {
                                              total_pop <- sum(df[[2]], na.rm = TRUE)
                                              if (total_pop == 0) return(NA)
                                              weighted_mean <- sum(df[[1]], na.rm = TRUE) / total_pop
                                              return(weighted_mean)
                                            },
                                            summarize_df = TRUE)
        icd_i_j_by_country <- exact_extract(rstack_foi_i_j_weighted, shp_country, fun = "sum")

        tp[[k]] <- list(country = tibble(ISO = shp_country@data$ISO,
                                         foi = foi_i_j_at_country,
                                         icd = icd_i_j_by_country,
                                         mod = fl_bkmod),
                        foiloc = tibble(loc_id = loc_foi@data$loc_id,
                                        pred_foi = foi_i_j_at_loc))
      }
      output <- list(
        country = tp[[1]]$country %>%
          left_join(tp[[2]]$country %>% rename(foi2 = foi, icd2 = icd, mod2 = mod), by = "ISO") %>%
          mutate(id_foi = j),
        foiloc = tp[[1]]$foiloc %>%
          left_join(tp[[2]]$foiloc %>% rename(pred_foi2 = pred_foi), by = "loc_id") %>%
          mutate(id_foi = j)
        )

      return(output)
    }

  print(proc.time() - ptm)
  stopCluster(cl)

  tp <- map_df(df_output, ~.x$country) %>%
    left_join(tibble(ISO = shp_country@data$ISO, suit = suit_i_at_country), by = "ISO")
  saveRDS(tp, file.path(wd_output, str_c("pred_incidence_country_suit", i, ".RDS")))
  message(str_c("Finished ", i, " out of", 101))
  tp2 <- map_df(df_output, ~.x$foiloc) %>%
    left_join(foi_est, by = c("id_foi" = "sim", "loc_id")) %>%
    left_join(suit_i_at_loc, by = "loc_id") %>%
    left_join(suit_i_at_loc2 %>% rename(suit2 = suit), by = "loc_id")
  saveRDS(tp, file.path(wd_output, str_c("pred_incidence_foiloc_suit", i, ".RDS")))
  
}
  
# Combine ----

tp <- map(c(1:9, 101), ~readRDS(file.path(wd_output, str_c("pred_incidence_country_suit", .x, ".RDS"))))
output <- bind_rows(tp, .id = "id_suit") %>% 
  left_join(tibble(ISO = shp_country@data$ISO, 
                   NAME = shp_country@data$NAME_ISO, 
                   pop = pop_at_country), 
            by = "ISO") %>%
  relocate(id_suit, id_foi, ISO, NAME, pop, suit, foi, icd, mod)

saveRDS(output, file = file.path(wd_output, "pred_incidence_country.RDS"))

