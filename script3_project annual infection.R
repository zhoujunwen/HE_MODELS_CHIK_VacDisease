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
select <- dplyr::select

# SERO-PREVALENCE 2025 & Cumulative incidence 2025-2050 ----

# > Proportion of people having incidence case in the year
ana <- readRDS(file = file.path(wd_data, "pred_incidence_country.RDS")) %>% 
  mutate(foi2 = icd / pop) %>% 
  select(id_suit, id_foi, ISO, propicd_annual = foi2)

# > Sampled incidence - base case
temp <- ana %>% 
  filter(id_suit != "101") %>% 
  group_by(ISO) %>% 
  sample_n(1000) %>% 
  ungroup %>% 
  arrange(ISO, id_suit, id_foi) %>% 
  group_by(ISO) %>% 
  mutate(id = row_number()) %>%
  ungroup %>% 
  select(id, ISO, propicd_annual)

# > Sampled incidence - sensitivity analysis
temp2 <- ana %>% 
  filter(id_suit == "101") %>% 
  select(ISO, propicd_annual) %>% 
  group_by(ISO) %>% 
  mutate(id = row_number() + 1000)  %>%
  ungroup

# > Different years of exposure in different countries
temp3 <- temp %>% distinct(ISO) %>% 
  mutate(firstyr_expo = case_when(
    # William M. de Souza 2024: https://www.sciencedirect.com/science/article/pii/S2667193X23002478
    # > America 2013
    ISO %in% c("BLZ", "BOL", "BRA", "COL", "CRI", "DOM", "ECU", "SLV", "GRD", "GTM",
               "GUY", "HTI", "HND", "MEX", "NIC", "PAN", "PRY", "PER", "LCA", "VEN") ~ 2013, 
    # Wahid 2017: https://pubmed.ncbi.nlm.nih.gov/28288924/
    ISO == "KHM" ~ 1961, # > Cambodia: KHM 1961 (The first case was identified in 1961)
    ISO %in% c("COD", "COG") ~ 2011, # > Democraci Repulic of Congo (COD), Republic of Congo (COG): Congo: The first outbreak occurred in 2011)
    ISO == "IND" ~ 1963, # > India: IND 1963 (The first outbreak occurred in Kolkata in 1963)
    ISO == "IDN" ~ 1972, # > Indonesia: IDN 1972 (In 1972, CHIKV was reported in East Sumatera, Kalimantan, Bali, Java, Sulawesi, and Flores)
    ISO == "MYS" ~ 1998, # > Malaysia: MYS 1998 (Several outbreaks have been reported in Port Klang (1998))
    ISO == "PHL" ~ 1965, # > Philippines: PHL 1965 (CHIKV was first identified in 1965)
    ISO == "SDN" ~ 1989, # > Sudan: SDN 1989 (About 12% of the population was positive for CHIKV in 1989)
    ISO == "THA" ~ 1960, # > Thailand: THA 1960 (46 000 suspected cases were reported in the 1960s)
    ISO %in% c("TCD", "DJI", "ETH") ~ 1952, # > Chad (TCD), Djibouti (DJI), Ethiopia (ETH) 1952; not reported, African country assume = 1952 in Zimbabwe earliest CHIK case in Africa
    TRUE ~ 0),
    firstyr_expo2 = ifelse(firstyr_expo == 2013, 2013, 1952)
  ) %>% 
  left_join(bind_rows(temp, temp2), by = c("ISO"))

# > cumulative seropos and annual prop of infection based on year of exposure

# > = Based on reported case 
tmp1 <- temp3 %>% 
  select(ISO, id, propicd_annual) %>% 
  expand_grid(yr_expo = 0:100) %>%
  arrange(ISO, id, yr_expo) %>% 
  group_by(id, ISO) %>% 
  mutate(propicd0 = 1-propicd_annual, 
         cumpropicd0 = cumprod(propicd0), # proportion of people not having infection at the end of the year
         seropos = 1 - lag(cumpropicd0, default = 1), #sero-postive at the beginning of the year
         propicd = propicd_annual * (1 - seropos) # 
  ) %>% 
  ungroup %>% 
  select(id, ISO, yr_expo, seropos, propicd)

# > Wrap-up
tmp2 <- temp3 %>% distinct(ISO, firstyr_expo) # Based on reported year is used

# > Projecting incidence ----

# Year of exposure since 
tmp3 <- map_df(unique(tmp2$firstyr_expo),
               ~expand_grid(tibble(firstyr_expo = .x, 
                                   yr_exp = .x:2050 - .x, 
                                   yr_act = .x:2050),
                            age = -100:100) %>% 
                 mutate(yr_age = age + yr_exp,
                        yr_exp2 = ifelse(age < 0, yr_age, yr_exp)) %>% 
                 filter(yr_age >= 0 & yr_age <100) %>% 
                 select(firstyr_expo, year = yr_act, age = yr_age, yr_expo = yr_exp2) %>%
                 arrange(year, age) %>% 
                 filter(year >= 2025))

# Population size over 2025-2050
temp <- read.csv(file.path(wd_data, "WPP2024_POP_F01_2_POPULATION_SINGLE_AGE_MALE - Medium variant.csv"))
temp2 <- read.csv(file.path(wd_data, "WPP2024_POP_F01_3_POPULATION_SINGLE_AGE_FEMALE - Medium variant.csv"))

# > Clean and summarize data by CHIK age group
tmp4 <- bind_rows(temp %>% mutate(sex = "m"),
                  temp2 %>% mutate(sex = "f")) %>%
  filter(Year %in% c(2025:2050) & ISO3.Alpha.code %in% unique(tmp1$ISO)) %>%
  select(country = ISO3.Alpha.code,
         year = Year,
         sex,
         starts_with("X")) %>%
  gather(key = "age", value = "n", -one_of(c("country", "year", "sex"))) %>%
  mutate(age2 = str_replace(str_remove(age, "X"), "\\.", "-"),
         age2 = ifelse(age2 == "100-", "100", age2),
         age3 = as.numeric(age2), 
         n2 = str_remove(trimws(n), " "),
         n3 = as.numeric(n2)) %>% 
  group_by(country, year, sex, age3) %>% 
  summarise(n = sum(n3), .groups = "drop") %>% 
  arrange(country, year, sex, age3) %>%
  select(country, year, sex, age = age3, n) %>% 
  filter(age < 100) %>% 
  mutate(n = ifelse(n == 0, 0.1, n * 1000)) # min > 0 is 0.5, set those 0 to 0.1

# Seroprevaelence over 2000-2024

# > One country by one country
tmp5 <- list()
temp0 <- unique(tmp2$ISO)
for(i in 1:length(temp0)){
  temp1 <- tmp1 %>% filter(ISO == temp0[i]) %>% select(id, yr_expo, seropos, propicd)
  temp2 <- tmp2 %>% filter(ISO == temp0[i]) %>% pull(firstyr_expo)
  temp3 <- tmp3 %>% filter(firstyr_expo == temp2) %>% select(-firstyr_expo)
  temp4 <- expand_grid(temp3, tibble(id = 1:1100))
  temp5 <- tmp4 %>% filter(country == temp0[i]) %>% select(year, sex, age, n) %>% spread(sex, n) 
  temp6 <- temp4 %>% left_join(temp1 , by = c("id", "yr_expo")) %>% 
    left_join(temp5, by = c("year", "age")) %>% 
    mutate(icd_f = f * propicd, icd_m = m * propicd,
           sero_f = f * seropos, sero_m = m * seropos)
  
  # Incidence 2025-2050
  temp6_0 <- temp6 %>% 
    mutate(ag = floor(age / 10)) %>% 
    group_by(id, year, ag) %>% 
    summarise_at(c("icd_f", "icd_m", "sero_f", "sero_m", "f", "m"), ~sum(.), .groups = "drop")
  temp6_1 <- temp6_0 %>% 
    filter(id %in% 1:1000) %>% 
    select(sim = id, year, sg_age = ag, f = icd_f, m = icd_m) %>%
    gather(key = "sex", value = "n_inf_bc", -one_of(c("sim", "sg_age", "year")))
  temp6_2 <- map_df(1:10,
                    ~temp6_0 %>% filter(id %in% 1001:1100) %>% mutate(id = id - .x * 100)) %>% 
    select(sim = id, year, sg_age = ag, f = icd_f, m = icd_m) %>%
    gather(key = "sex", value = "n_inf_sa", -one_of(c("sim", "sg_age", "year")))
  tmp5[[i]] <- left_join(temp6_1, temp6_2, by = c("sim", "year", "sex", "sg_age"))  %>%
    mutate(country = temp0[i])
  print(i)
}

tp <- bind_rows(tmp5) %>% 
  relocate(sim, country, year, sex, sg_age, n_inf_bc, n_inf_sa) %>% 
  arrange(sim, country, year, sex, sg_age) %>% 
  ungroup()
saveRDS(tp, file.path(wd_output, "ana_epi_inf_by_age_sex.rds"))
