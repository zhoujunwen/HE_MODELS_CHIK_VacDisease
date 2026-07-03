# HE_MODELS_CHIK_VacDisease

Health-economic burden of chikungunya infection and potential cost-effectiveness of preventive vaccination in 31 countries

This folder contains the example codes (and some with small set of data) used in a vaccine impact modelling study assessing the health economic benefit of vaccination. 

Overall, the model consists of six components. 
1. Global CHIKV vectors suitability was estimated using both random forest and boosted regression trees (BRT), vector occurrence data and environmental covariates. 
2. Global CHIKV suitability was estimated using BRT, recorded disease occurrence data and environmental covariates. 
3. Historical patterns of the force of infection (FoI) in different areas of the world estimated from seroprevalence data were taken from Kang et al (https://pubmed.ncbi.nlm.nih.gov/38342105/). 
4. The relationship between suitability and FoI was estimated using a logistic growth model. 
5. This relationship and estimates of the population at risk, derived from antibody prevalence estimates and a temperature-based mask to exclude areas deemed unlikely suitable for sustained mosquito survival and virus replication, were used to project average annual CHIKV infection incidence forward in time from 2025 to 2050. 
6. A decision-analytic model was used to estimate country-specific health-economic outcomes of CHIKV infection with and without CHIKV vaccine administration.

This folder contains four scripts corresponding to the above analyses
1. Estimate CHIKV suitability
2. Predict CHIKV FOI and incidence 
3. Project annual infection
4. Estimate health-economic burden and vaccine impacts

Relevant publication: 
Zhou J, Salant N, Radoykova HS, Messina J, Wint W, Longbottom J, Holohan KM, Dinkel KA, Sytsma MLT, Torkelson AA, Cavalcanti LPG, Hollingsworth TD, Lord J, Smith DRM, Pouwels K. 
Health-economic burden of chikungunya infection and potential cost-effectiveness of preventive vaccination in 31 countries. 2026. Medrxiv link. 
