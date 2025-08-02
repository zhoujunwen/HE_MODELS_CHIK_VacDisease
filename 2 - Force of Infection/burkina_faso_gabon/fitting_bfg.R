library('ggplot2')
library('purrr')
library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library('tidyverse')
library('reshape2')
library('bayesplot')
library('loo')
library('pracma')
library('cowplot')
library('grid')
library('gridExtra')
library('Hmisc')
library('dplyr')
# library(vscDebugger) # uncomment and run in VSCode
library('epitrix')
library('gsubfn') # obtain limits from applying cut function


source('utils.R')

# this file is created using preprocessing.R which accesses `open_data`
all_df = get(load('burkina_faso_gabon/data/df_burkina_faso_gabon_2015.RData'))

# Stan model to be compiled before running 
# MContNormal    <- stan_model('mod/target_bimod_prior.stan') 
# bimodal prior with foi (working sims)
# MContNormal    <- stan_model('burkina_faso_gabon/target_bimod_prior_forloop.stan')
# serofoi model 
MContNormal    <- stan_model('burkina_faso_gabon/tv_normal_log.stan')
cat('Model compiled successfully. \n')
#load('MContNormal.RDS') # will work if previously compiled

# regions to be iterated over
all_regions = all_df$Region %>% unique

#### Parameters 
niters =  30e3#4e3 # 30e3   # massively reduce niters 
nwarmup =  20e3#2e3 # 20e3  # 5k warmup + 10k iters -- or 2k and 4k  
delta_adpt = 0.92 # 0.92 for convergence issues suggested by src code 
nchains = 4

# make sure to add '/' at the end 
dest_dir = 'burkina_faso_gabon/res_bfg/'

fit_res_errs <- map(all_regions, function(region_i){
    t0 <- Sys.time()
    # subset data
    dat <- filter(all_df, Region == region_i) %>% 
        arrange(., age_mean_f) %>% 
        mutate(., birth_year = Year - age_mean_f)
    # generate code for saved files
    survey_code = paste(dat$Country[1], dat$Region[1], dat$Year[1], sep='_')
    cat(paste("Fitting to: ______ ______ ", survey_code,"\n", sep=""))
    
    # exposure matrix 
    # for some reason it is needed in the global environment 
    yexpo <<- make_yexpo(dat)
    # actual years 
    RealYexpo <- (min(dat$birth_year):dat$Year[1])[-1]
    # create expo matrix = specify the FOIs that every age group experiences
    ExposureMatrix <- get_exposure_matrix(dat, yexpo)
    
    # data fet to stan model sampler
    stan_data <- list(
        Nobs   = nrow(dat),
        Npos   = dat$N_seropos,
        Ntotal = dat$N_tot,
        Age    = dat$age_mean_f,
        Ymax   = max(yexpo),
        AgeExpoMatrix = ExposureMatrix,
        NDecades = 0
        )
    # sampler fits model to data 
    fit_tmp <- safe_quiet_fit(MContNormal,
        data = stan_data,
        iter = niters,
        chains = nchains,
        init = fInit, # feed function to initialise the system
        warmup = nwarmup,
        verbose = FALSE,
        refresh = 0,
        control = list(
            adapt_delta = delta_adpt, # default 0.8
            # made 0.9 or 0.92 for poorer convergence
            max_treedepth = 12
        ), # cannot find default in docs
        # internet val is 10
        seed = "12345",
        thin = 1
    ) %>% tidy_up
    fit <- fit_tmp$result
    if (!is.null(fit_tmp$output)) {
        if (!is.null(fit_tmp$error)) {
            fit_tmp$error <- c(fit_tmp$error, fit_tmp$output)
        } else {
            fit_tmp$error <- fit_tmp$output
        }
    }
    # examining parameters if errors are thrown out 
    #pairs(fit) # errors due to large num of params 
    # need to choose which params to visualise
    if(!is.null(fit@sim$samples)) { # check if fitting was successful 
        # leave-one-out cross-validation (approximate)
        loo_fit_tmp <- safe_quiet_loo(fit, save_psis = TRUE, 'logLikelihood') %>% tidy_up
        loo_fit <- loo_fit_tmp$result
        if (!is.null(loo_fit_tmp$warnings)) {
            fit_tmp$loo_warnings <- loo_fit_tmp$warnings
        } else { 
            fit_tmp$loo_warnings = NULL 
            }
        if (!is.null(loo_fit_tmp$loo_messages)) {
            fit_tmp$loo_messages <- loo_fit_tmp$messages
        } else { 
            fit_tmp$loo_messages = NULL 
            }
        if (!is.null(loo_fit_tmp$errors)) {
            fit_tmp$errors <- loo_fit_tmp$errors
        } else { 
            fit_tmp$errors = NULL 
            }
        # Extract samples from a fitted model
        foi <- rstan::extract(fit, 'foi', inc_warmup = FALSE)[[1]]
        # measures of central tendency + CI
        foi_cent_est <- data.frame(year  = RealYexpo,
                                lower = apply(foi, 2, function(x) quantile(x, 0.1)),
                                upper = apply(foi, 2, function(x) quantile(x, 0.9)),
                                medianv = apply(foi, 2, function(x) quantile(x, 0.5)))
        # sample from the FOI -- see supplementary material
        # used for the burden model, can ignore for CHIK purposes 
        # https://royalsocietypublishing.org/doi/suppl/10.1098/rstb.2022.0278
        foi_post_1000s <- dplyr::sample_n(as.data.frame(foi), size = 1000)
        colnames(foi_post_1000s) <- RealYexpo
        
        ## structure of res file 
        res <- list(fit=fit,
                    stan_data = stan_data,
                    RealYexpo = RealYexpo,
                    yexpo     = yexpo,
                    n_iters   = niters,
                    n_thin    = 1,
                    n_warmup  = nwarmup,
                    model     = 'time-var-FOI',
                    delta     = delta_adpt,
                    mtreed    = 12,
                    loo_fit   = loo_fit,
                    foi_cent_est   = foi_cent_est,
                    foi_post_1000s  = foi_post_1000s)

        # save results as Rdata files 
        title_data = paste0(dest_dir, survey_code, ".RData")
        save(res, file=title_data)
        
        # determine max lambda for plotting 
        foi_mod <- rstan::extract(res$fit, 'foi', inc_warmup = FALSE)[[1]]
        max_lambda <-  (as.numeric(quantile(foi_mod, 0.95))) * 1.3
        
        # combined plots object 
        pcc_tmp <- safe_quiet_fCombinedPlots(
            res = res, dat = dat, lambda_sim = NA,
            max_lambda = max_lambda, size_text = 12
        ) %>% tidy_up
        pcc <- pcc_tmp$result
        fit_tmp$plot_warnings <- pcc_tmp$warnings
        fit_tmp$plot_messages <- pcc_tmp$messages
        if (class(pcc)!='logical'){
            plot_arrange <- grid.arrange(pcc$plots[[1]], 
                                        pcc$plots[[2]], 
                                        pcc$plots[[3]], 
                                        pcc$plots[[4]], 
                                        nrow = 4, heights = c(1.5, 1, 1, 1))
            
            title = paste(dest_dir, survey_code, ".png", sep='')
            ggsave(title, plot=plot_arrange, dpi=330, width=5, height=8, units='in')
        }
        tf <- Sys.time()
        time_taken <- tf-t0
        print(time_taken)
        
    } else { # if no fit@sim$samples, return no model 
        cat('ERROR:  NO MODEL FIT FOR  ')
        cat(survey_code)
        cat('\n')
        res = NA
        }

    fit_tmp$region = region_i 
    return(fit_tmp)
})

# save logs in desired format 
write_log_json(fit_res_errs, dest_dir=dest_dir, save_csv=FALSE) # can additionally save CSV
# write_log_txt(fit_res_errs, dest_dir=dest_dir)


# this produces some ggplot warnings, 
# likely from plot differences that i have set up 
# e.g. the more narrow time frame 
fit_res_errs

# check errors and warnings like so  
#### alternatively check saved log file 
map(fit_res_errs, c('error')) %>% walk(cat) # %>% compact uncomment to discard NULLs
map(fit_res_errs, c('warnings')) 
map(fit_res_errs, c('messages')) 
map(fit_res_errs, c('loo_warnings')) %>% compact
map(fit_res_errs, c('plot_warnings')) %>% compact
map(fit_res_errs, c('plot_messages')) %>% compact 
    