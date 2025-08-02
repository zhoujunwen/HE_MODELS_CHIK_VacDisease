
# Start ----

## Load data ----

# Clear Environment
rm(list=ls())

# Set work directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # the directory of this script

# Include package
library("tidyverse")
select<- dplyr::select

# Load data
load("ana4_data.Rdata")

## Run model ----

# > Scaling the infection number over years adjusting for seroprevalence and population size

vec_country <- sort(ana_sero$country)
vec_years <- 2025:2050

tmp <- rep(list(), length(vec_country))
for(i in 1:length(vec_country)){
  country_i <- vec_country[i]
  temp <- ana_inf %>% filter(country == country_i) 
  temp2 <- ana_sero %>% filter(country == country_i) %>% pull(weighted_avg)
  temp3 <- ana_pop %>% filter(country == country_i)
  vec_N = temp3$nt
  vec_B = temp3$nb
  vec_D = temp3$nd
  tmp_i <- list()
  
  N_0 = vec_N[1]
  P_0 = temp2
  if(length(P_0) > 1){stop("multiple baseline seroprevalence estimates")}
  S_0 = P_0 * N_0
  for(j in 1:1000){
    
    ### Baseline annual values: population size (N), seroprevalence (P), number seropositive (S), incidence (Inc)
    Inc_0 = temp$incidence[j]
    
    ### empty vectors
    vec_P = c()
    vec_S = c()
    vec_Inc = c()
    
    ### first year
    N_t = N_0
    P_t = P_0
    S_t = S_0
    Inc_t = Inc_0 * (1 - P_0)
    
    qounter_year = 0
    for(year_j in vec_years){
      
      qounter_year = qounter_year + 1
      
      # current year's values
      vec_S[qounter_year] = S_t
      vec_P[qounter_year] = P_t
      vec_Inc[qounter_year] = Inc_t
      
      # this year's demography
      N_current = vec_N[qounter_year]
      B_t = vec_B[qounter_year]
      D_t = vec_D[qounter_year]
      
      # Scaling both population size and seroprevalence
      
      # > POP SIZE UPDATE
      N_t = N_current + B_t - D_t
      
      # > S AND P UPDATE
      S_t = S_t + Inc_t - P_t * D_t
      P_t = S_t / N_t
      
      # incidence next year: scale by population size and seroprevalence
      Inc_t = Inc_0 * (N_t / N_0) * (1 - P_t)
    }
    tmp_i[[j]] <- tibble(sim = j, 
                            country = country_i,
                            year = 2025:2050,
                            incidence = vec_Inc)
  }
  tmp[[i]] <- bind_rows(tmp_i)
  print(i)
}

output <- bind_rows(tmp) %>% rename(n_inf = incidence)
name_out <- "ana_epi_inf.rds" 
saveRDS(output, name_out)
