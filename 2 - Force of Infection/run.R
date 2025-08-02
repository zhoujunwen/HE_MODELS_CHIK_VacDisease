

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd()
library("ggplot2")
library("purrr")
library("rstan")
options(mc.cores = parallel::detectCores()) # this was suggested
rstan_options(auto_write = TRUE)            # by rstan docs
library("tidyverse")
library("reshape2")
library("bayesplot")
library("loo")
library("pracma")
library("cowplot")
library("grid")
library("gridExtra")
library("Hmisc")
library("dplyr")
# library(vscDebugger) # uncomment only in VSCode
library("epitrix")
library("gsubfn") # obtain limits from applying cut function

# fitting functions and error handling helpers
source("utils.R")
# function with for-loop to fit over all regions
source("fitting.R")


# this large data file is created using preprocessing.R which accesses `open_data`
# 'burkina_faso_gabon/data/df_burkina_faso_gabon_2015.RData'
data_dir <- "data/all_data_kang.RData"
# data_dir <- "data/all_data_open_access.RData"
df_all_regions <- get(load(data_dir))

# results directory
dest_dir <- "res_kang" #' res_debug/'# 'res/'
if (!dir.exists(dest_dir)) {
    dir.create(dest_dir)
}

# Stan model to be compiled before running 
# time var with bimodal prior   
# MContNormal    <- stan_model('mod/target_bimod_prior_forloop.stan')
# serofoi fast epidemic model
### log-scale FOI with normal prior and cauchy st dev 
MContNormal    <- stan_model('mod/tv_normal_log.stan')
### DEBUG certain fitting defects using a stronger prior 
# stronger prior = use Cauchy range [0, 0.4] 
# MContNormal    <- stan_model('mod/tv_normal_log_cauchy_04.stan')
cat('Model compiled successfully. \n')
#load('MContNormal.RDS') # will work if previously compiled

# regions to be iterated over
# have to be unique 
all_regions <- df_all_regions$study_id %>% unique()
# debugging 
# all_regions = 'Kumbo' # 'Guadeloupe' # 'Kumbo'


#### Parameters 
# niters <-     30e3 # 4e3 # massively reduce niters for rough results 
# nwarmup <-    20e3 # 2e3 # e.g. 5k warmup + 10k iters -- or 2k and 4k as shown 
# delta_adpt <- 0.92 # 0.92 but 0.95 for convergence issues suggested by src code, up to 0.98
# nchains <-    4 # when there are "errors in some chains" change to 1 to debug

niters <-     4e3 # massively reduce niters for rough results 
nwarmup <-    2e3 # e.g. 5k warmup + 10k iters -- or 2k and 4k as shown 
delta_adpt <- 0.92 # 0.92 but 0.95 for convergence issues suggested by src code, up to 0.98
nchains <-    4 # when there are "errors in some chains" change to 1 to debug


# if n_chains > 1, opens a viewer window in RStudio 
# or a browser tab in vscode
# but neither work for my setup 
# see https://discourse.mc-stan.org/t/rstudio-viewer-with-rstan-running-in-parrallel-with-rstudio-1-3/15615
fit_res_errs <- fit_to_regions(
    all_regions,    # regions or region IDs to be looped over 
    save_logs=T     # "dest_dir/log_YYMMDD_HHMMSS.json" (optionally csv, see utils)
    )
# NOTE code produces some ggplot warnings, 
#### likely from plot differences that i have set up 
#### e.g. the more narrow time frame plotted 
# but they can be ignored 


# check model fitting errors and warnings like so  
#### alternatively check saved log file 
map(fit_res_errs, c('error')) %>% walk(cat) # %>% compact # discard NULLs
map(fit_res_errs, c('warnings')) 
map(fit_res_errs, c('messages')) 
map(fit_res_errs, c('loo_warnings')) %>% compact
map(fit_res_errs, c('plot_warnings')) %>% compact
map(fit_res_errs, c('plot_messages')) %>% compact 
    


# PPC1    <- fCombinedPlots(res=mod_1, dat, lambda_sim, max_lambda)
# mod_1$prev_expanded <- PPC1$prev_expanded
# parrange_1 <- vertical_plot_arrange_per_model(PPC1)
# 
# pcc$prev_expanded
# 
# parrange = vertical_plot_arrange_per_model(pcc)

# SAVE
# save(fit_res_errs, file='res/fit_res_errs.RData')
# res_saved = get(load('res/fit_res_errs.RData')) # it's a large file, might fail

