
# Start ----

# Set work directory (Where you load the data)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # the directory of this script

# Include package
library("tidyverse")
select<- dplyr::select

# Prepare data ----

# Load data
load("ana5_data.RData")

# Prepare scenario
ana_input <- tibble(symp = "bc", dw_chronic = "bc", scn = 1)

# Prepare function ----

# Run simulation
f_run_simulation <- function(input = ana_input, opt_out = TRUE, opt_diff = FALSE){
  
  # input <- ana_input
  mod_epi_oth <- ana_epi_oth %>% 
    rename(p_symp = str_c("p_symp_", input$symp), 
           dw_chronic = str_c("dw_chronic_", input$dw_chronic)) %>% 
    select(-ends_with("_bc"), -ends_with("_sa"), -ends_with("_sa1"), -ends_with("_sa1"))
  mod_scn <- ana_scn[[input$scn]]
  
  vac_i <- mod_scn$vac
  vac_eff_i <- mod_scn$eff
  
  # All parameters
  tp <- ana_epi_inf %>%
    left_join(mod_epi_oth, by = c("sim", "sg_age", "sex")) %>% 
    left_join(ana_epi_ly, by = c("country", "year", "sg_age", "sex")) %>% 
    left_join(ana_eco_country_year_sg, by = c("country", "year", "sg_age", "sex")) %>% 
    left_join(ana_eco_country, by = "country") %>% 
    left_join(vac_i, by = c("country", "year", "sg_age", "sex")) 
  
  # N
  n_inf <- tp$n_inf
  n_symp <- with(tp, n_inf * p_symp * (prop_vac_scn * vac_eff_i + (1 - prop_vac_scn))) 
  n_severe <- n_symp * tp$p_severe
  n_mild <- n_symp - n_severe
  n_death <- n_severe * tp$p_death
  n_chronic <- (n_symp - n_death) * tp$p_chronic
  
  # DUR
  y_mild <-  n_mild * tp$dur_mild
  y_prehosp <- n_severe * tp$dur_severe_prehosp
  y_hosp <- n_severe * tp$dur_hosp
  y_chronic <- n_chronic * tp$dur_chronic
  
  # DALY by state
  daly_mild <- y_mild * tp$dw_mild
  daly_severe <- y_prehosp * tp$dw_severe_prehosp + y_hosp * tp$dw_hosp
  daly_chronic <- y_chronic * tp$dw_chronic
  daly_death <- n_death * tp$ly
  daly_death_disc <- n_death * ((1/0.035)*(1-exp(-0.035*(tp$ly)))) 
  
  # DALY total
  daly_tot <- daly_mild + daly_severe + daly_chronic + daly_death
  daly_tot_disc <- daly_tot - daly_death + daly_death_disc # only for costing
  
  # Care
  n_op_mild <- n_mild * tp$prop_care_mild
  n_op_chronic <- n_chronic * tp$prop_care_chronic
  n_op_chronic_mild_op0 <- (n_mild - n_op_mild) * tp$prop_care_chronic

  # Economics
  
  # > Costs
  cost_care_mild <- n_op_mild * tp$cost_care_mild 
  cost_care_severe <- n_severe * tp$cost_care_severe
  cost_care_chronic <- n_op_chronic * tp$cost_care_chronic 
  cost_care <- cost_care_mild + cost_care_severe + cost_care_chronic
  
  # > Costs OOP
  cost_care_oop_mild <- n_op_mild * tp$cost_care_mild_oop
  cost_care_oop_severe <- n_severe * tp$cost_care_severe_oop
  cost_care_oop_chronic <- n_op_chronic * tp$cost_care_chronic_oop
  cost_care_oop <- cost_care_oop_mild + cost_care_oop_severe + cost_care_oop_chronic
  
  # > CHE & IHE
  n_che_mild <- n_op_mild * tp$prop_che_mild
  n_che_severe <- n_severe * tp$prop_che_hosp
  n_che_chronic <- n_op_chronic_mild_op0 * tp$prop_che_chronic
  n_che = n_che_mild + n_che_chronic + n_che_severe
  n_che_work = n_che * tp$prop_working
  
  n_ihe_mild <- n_op_mild * tp$prop_ihe_mild 
  n_ihe_severe <- n_severe * tp$prop_ihe_hosp
  n_ihe_chronic <- n_op_chronic_mild_op0 * tp$prop_ihe_chronic
  n_ihe = n_ihe_mild + n_ihe_chronic + n_ihe_severe
  n_ihe_work = n_ihe * tp$prop_working
  
  # > Costs DALY
  cost_daly <- daly_tot * tp$cost_daly
  cost_daly_disc <- daly_tot_disc * tp$cost_daly
  cost_daly_mild <- daly_mild * tp$cost_daly
  cost_daly_severe <- daly_severe * tp$cost_daly
  cost_daly_chronic <- daly_chronic * tp$cost_daly
  cost_daly_death_disc <- daly_death_disc * tp$cost_daly
  
  # > Costs VSL
  cost_vsl <- n_death * tp$cost_vsl
  cost_vsly <- n_death * tp$cost_vsl / tp$ly
  
  # > PROD LOSS
  y_prodloss_mild <- y_mild * tp$prop_working * ana_eco$prodloss$prod_impact$mild
  y_prodloss_severe <- (y_prehosp + y_hosp) * tp$prop_working * ana_eco$prodloss$prod_impact$severe
  y_prodloss_chronic <- y_chronic * tp$prop_working * ana_eco$prodloss$prod_impact$chronic
  y_prodloss_death <- n_death * tp$ly_working
  y_prodloss_tot <- y_prodloss_mild + y_prodloss_severe + y_prodloss_chronic + y_prodloss_death
  y_prodloss_death_disc <- n_death * tp$ly_working_disc
  
  cost_prodloss <- y_prodloss_tot * tp$cost_prod
  cost_prodloss_disc <- (y_prodloss_mild + y_prodloss_severe + y_prodloss_chronic + y_prodloss_death_disc) * tp$cost_prod
  
  cost_prodloss_mild <- y_prodloss_mild * tp$cost_prod
  cost_prodloss_severe <- y_prodloss_severe * tp$cost_prod
  cost_prodloss_chronic <- y_prodloss_chronic * tp$cost_prod
  cost_prodloss_death_disc <- y_prodloss_death_disc * tp$cost_prod
  
  # > Societal costs
  cost_societal <- cost_care + cost_daly + cost_prodloss
  cost_societal_vsly <- cost_care + cost_vsly + cost_prodloss
  cost_societal_disc <- cost_care + cost_daly_disc + cost_prodloss_disc
  cost_societal_vsly_disc <- cost_care + cost_vsly + cost_prodloss_disc
  
  tmp <- tp %>% 
    select(sim, country, year, sex, sg_age) %>% 
    bind_cols(
      tibble(n_inf = n_inf, 
             n_symp = n_symp,
             n_severe = n_severe,
             n_death = n_death,
             n_chronic = n_chronic,
             daly_mild = daly_mild,
             daly_severe = daly_severe,
             daly_chronic = daly_chronic,
             daly_death = daly_death,
             daly_tot = daly_tot,
             n_che = n_che,
             n_che_work = n_che_work,
             n_ihe = n_ihe,
             n_ihe_work = n_ihe_work,
             cost_care = cost_care,
             cost_prodloss = cost_prodloss,
             cost_vsl = cost_vsl,
             cost_vsly = cost_vsly,
             cost_daly = cost_daly,
             cost_societal = cost_societal,
             cost_societal_vsly = cost_societal_vsly,
             cost_care_disc = cost_care,
             cost_prodloss_disc = cost_prodloss_disc,
             cost_vsl_disc = cost_vsl,
             cost_vsly_disc = cost_vsly,
             cost_daly_disc = cost_daly_disc,
             cost_societal_disc = cost_societal_disc,
             cost_societal_vsly_disc = cost_societal_vsly_disc, 
             cost_care_mild_disc = cost_care_mild,
             cost_care_severe_disc = cost_care_severe,
             cost_care_chronic_disc = cost_care_chronic,
             cost_prodloss_mild_disc = cost_prodloss_mild,
             cost_prodloss_severe_disc = cost_prodloss_severe,
             cost_prodloss_chronic_disc = cost_prodloss_chronic,
             cost_prodloss_death_disc = cost_prodloss_death_disc,
             cost_daly_mild_disc = cost_daly_mild,
             cost_daly_severe_disc = cost_daly_severe,
             cost_daly_chronic_disc = cost_daly_chronic,
             cost_daly_death_disc = cost_daly_death_disc)
      )%>% 
    mutate_at(str_c("cost_",
                    c("care", "prodloss", "daly", "vsl", "vsly", "societal", "societal_vsly",
                      "care_mild", "care_severe", "care_chronic", 
                      "prodloss_mild", "prodloss_severe", "prodloss_chronic", "prodloss_death",
                      "daly_mild", "daly_severe", "daly_chronic", "daly_death"),
                    "_disc"), 
              ~. / (1+0.035)^(year - 2025))
  tmp2_1 <- tmp %>% filter(sex == "f")
  tmp2_2 <- tmp %>% filter(sex == "m")
  output <- bind_cols(tmp2_1 %>% select(sim, country, year, sg_age),
                      map2(tmp2_1[,-(1:5)], tmp2_2[,-(1:5)], ~.x + .y))
  
  # Burden in the scenario
  if(opt_out == TRUE){
    saveRDS(output, 
            str_c("rst_symp", input$symp, 
                  "_dw", input$dw_chronic, 
                  "_scn", input$scn, "_diff0.rds"))
  }
  
  # Difference compared to the base-case scenario
  if(opt_diff == TRUE){ 
    output_base <- readRDS(
      str_c("rst_symp", input$symp, 
            "_dw", input$dw_chronic, 
            "_scn", 1, "_diff0.rds"))
    output <- bind_cols(
      output %>% select(sim, country, year, sg_age),
      map2(output[,-(1:4)], output_base[,-(1:4)], ~.x - .y))
    saveRDS(output, 
            str_c("rst_symp", input$symp, 
                  "_dw", input$dw_chronic, 
                  "_scn", input$scn, "_diff1.rds"))
    }
}

# Summarize the results from the simulation
f_sum_simulation <- function(input = ana_input, opt_diff = 0, opt_global = TRUE, opt_country = FALSE){

  ana <- readRDS(str_c("rst_symp",input$symp, 
                       "_dw", input$dw_chronic, 
                       "_scn", input$scn, 
                       "_diff", opt_diff, ".rds"))
  
  # Summary across all the countries
  if(opt_global == TRUE){
    temp <- ana %>%
      select(-country) %>%
      group_by(sim, year, sg_age) %>%
      summarise_all(~sum(.)) %>%
      ungroup
    tmp1_1 <- temp
    tmp1_2 <- temp %>%
      select(-year) %>%
      group_by(sim, sg_age) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(year = 2024)
    tmp1_3 <- temp %>%
      select(-sg_age) %>%
      group_by(sim, year) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(sg_age = 0)
    tmp1_4 <- tmp1_3 %>%
      select(-year) %>%
      group_by(sim) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(year = 2024, sg_age = 0)
    tmp2 <- bind_rows(tmp1_1, tmp1_2, tmp1_3, tmp1_4) %>%
      select(-sim) %>%
      group_by(year, sg_age) %>%
      summarise_all(list(m = mean,
                         l = ~quantile(., probs = 0.025),
                         h = ~quantile(., probs = 0.975))) %>%
      ungroup %>%
      arrange(year, sg_age)
    output <- map_df(c("m", "l", "h"),
                     ~tmp2 %>%
                       select(c("year", "sg_age"), ends_with(str_c("_", .x))) %>%
                       gather(key = "est", value = "val", -one_of(c("year", "sg_age"))) %>%
                       mutate(out = str_sub(est, 1, -3),
                              para = .x) %>%
                       select(-est)) %>%
      spread(key = "para", value = "val")
    name_out <- str_c("sum_symp",input$symp, 
                      "_dw", input$dw_chronic, 
                      "_scn", input$scn, 
                      "_diff", opt_diff, 
                      "_by_all.rds")
    saveRDS(output, name_out)
  }
  
  # Summary at the country-level
  if(opt_country == TRUE){
    tmp1_1 <- ana
    tmp1_2 <- ana %>%
      select(-year) %>%
      group_by(sim, country, sg_age) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(year = 2024)
    tmp1_3 <- ana %>%
      select(-sg_age) %>%
      group_by(sim, country, year) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(sg_age = 0)
    tmp1_4 <- tmp1_3 %>%
      select(-year) %>%
      group_by(sim, country) %>%
      summarise_all(~sum(.)) %>%
      ungroup %>%
      mutate(year = 2024, sg_age = 0)
    tmp2 <- bind_rows(tmp1_1, tmp1_2, tmp1_3, tmp1_4) %>%
      select(-sim) %>%
      group_by(country, year, sg_age) %>%
      summarise_all(list(m = mean,
                         l = ~quantile(., probs = 0.025),
                         h = ~quantile(., probs = 0.975))) %>%
      ungroup %>%
      arrange(country, year, sg_age)
      output <- map_df(c("m", "l", "h"),
                       ~tmp2 %>%
                         select(c("country", "year", "sg_age"), ends_with(str_c("_", .x))) %>%
                         gather(key = "est", value = "val", -one_of(c("country", "year", "sg_age"))) %>%
                         mutate(out = str_sub(est, 1, -3),
                                para = .x) %>%
                         select(-est)) %>%
        spread(key = "para", value = "val")
      name_out <- str_c("sum_symp", input$symp, 
                        "_dw", input$dw_chronic, 
                        "_scn", input$scn, 
                        "_diff", opt_diff, 
                        "_by_country.rds")
      saveRDS(output, name_out)
  } 
}

# Run simulation ----

# Health-economic burden across countries and at each country 
ana_input <- tibble(symp = "bc", dw_chronic = "bc", scn = 1) 
f_run_simulation(input = ana_input, opt_out = TRUE, opt_diff = FALSE)
f_sum_simulation(input = ana_input, opt_diff = 0, opt_global = TRUE, opt_country = TRUE)

# Vaccination impact across countries: 3 vaccine efficacy x 7 vaccination compaign
for(i in 2:22){
  ana_input <- tibble(symp = "bc", dw_chronic = "bc", scn = i) 
  f_run_simulation(input = ana_input, opt_out = FALSE, opt_diff = TRUE)
  f_sum_simulation(input = ana_input, opt_diff = 1, opt_global = TRUE, opt_country = FALSE)
}

# Vaccination impact at each country: vaccine 70% efficacy using population wide / routine / combined complain
scn_interest <- c(2,4,6) + 1 + 7 
# 2: middle of population wide compaign
# 4: routine compaign
# 6: middle of combine compaign
# + 1 -> shift away from baseline burden
# + 7 -> shift away from 90% efficacy
for(i in 1:3){
  ana_input <- tibble(symp = "bc", dw_chronic = "bc", scn = scn_interest[i]) 
  f_sum_simulation(input = ana_input, opt_diff = 1,opt_global = FALSE, opt_country = TRUE)
  print(i)
}

