# HE_MODEL_CHIK_DiseaseVac

Health and economic burden of chikungunya infection and potential benefits of vaccination in 32 countries: a vaccine impact modelling study

This folder contains the codes and small set of data used in a vaccine impact modelling study assessing the health economic benefit of vaccination. 

Overall, the model consists of five main components. 
1. Global CHIKV suitability was estimated using boosted regression trees (BRT), recorded disease occurrence data and environmental covariates. 
2. Historical patterns of the force of infection (FoI) in different areas of the world were estimated from seroprevalence data. 
3. The relationship between suitability and FoI was estimated using a logistic growth model. 
4. This relationship and estimates of the population at risk, derived from antibody prevalence estimates and a temperature-based mask to exclude areas deemed unlikely suitable for sustained mosquito survival and virus replication, were used to project average annual CHIKV infection incidence forward in time from 2025 to 2050. 
5. A decision-analytic model was used to estimate country-specific health-economic outcomes of CHIKV infection with and without CHIKV vaccine administration.

Therefore, this folder contains five folders corresponding to the analyses used in each component
1. Estimate CHIKV suitability
2. Estimate Force of infection (Forked from https://github.com/sedaradoykova/)
3. Map suitability to FOI
4. Project annual infection
5. Estimate healt-economic burden and vaccine impacts

Relevant publication: 
Zhou J, Salant N, Radoykova HS, Messina J, Wint W, Longbottom J, Holohan KM, Sytsma MLT, Torkelson AA, Cavalcanti LPG, Mukandavire C, Lord J, Smith DRM, Hollingsworth TD, Pouwels K. 
Health and economic burden of chikungunya infection and potential benefits of vaccination in 32 countries: a vaccine impact modelling study. 2025. 
Medrxiv link. 
