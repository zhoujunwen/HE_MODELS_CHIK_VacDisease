library("jsonlite")

########################
#### MODEL HELPERS ####
#######################

get_exposure_matrix <- function(dat, yexpo) {
    age_class <- dat$age_mean_f
    ly <- length(yexpo)
    # exposure mat has dims n_age_groups X n_years
    exposure <- matrix(0, nrow = length(age_class), ncol = ly)
    # fill out such that oldest get 1 across years (i.e. exposed to all FOIs)
    # youngest only get exposed to most recent FOIs while they are alive
    for (k in 1:length(age_class)) {
        ind_min <- ly - age_class[k] + 1
        # ifelse is important for negative indexing error
        ind_min <- ifelse(ind_min < 0, 1, ind_min)
        exposure[k, ind_min:ly] <- 1
    }
    exposure_output <- exposure
    return(exposure_output)
}

make_yexpo <- function(dat) {
    yexpo <- (seq_along(min(dat$birth_year):dat$Year[1]))
    yexpo <- yexpo[-length(yexpo)]
}



# function initialises all parameters
#### FOIs = 0.01 but should start chains from different starting points (?)
# fInit <- function() {
#     list(foi=rep(0.01, length(yexpo)))
# } 
### for target_bimod_prior_forloop.stan
# this is for log-scale FOI in the explosive epidemic model 
fInit <- function() {
    list(log_foi = rep(-3, length(yexpo)))
}


##########################
#### PLOTTING HELPERS ####
##########################
## these were written by 
## Ledien et al 

vertical_plot_arrange_per_model <- function(PPC) {
    pp <- grid.arrange(PPC$plots$plot_data,
        PPC$plots$plot_prev,
        PPC$plots$plot_foi,
        PPC$plots$plot_rhats,
        nrow = 4,
        heights = c(1.5, 1, 1, 1)
    )
    return(pp)
}

extract_summary_mod <- function(res, dat) {
    model_name <- res$model
    #------- Loo estimates

    loo_fit <- res$loo_fit
    if (sum(is.na(loo_fit)) < 1) {
        lll <- as.numeric((round(loo_fit$estimates[1, ], 2)))
    } else {
        lll <- c(-1e10, 0)
    }

    # browser() # Okay! Aquí obtengo bien lll con elpd

    summary_mod <- data.frame(
        model = res$model,
        country = dat$Country[1],
        region = dat$Region[1],
        year = dat$Year[1],
        n_sample = sum(dat$N_tot),
        n_agec = length(dat$age_mean_f),
        n_iter = res$n_iters,
        # performance = "_____",
        elpd = lll[1],
        se = lll[2],
        converged = NA
    )
    rhats <- get_table_rhats(res)
    if (any(rhats$rhat > 1.1) == FALSE) {
        summary_mod$converged <- "Yes"
    }

    return(summary_mod)
}



get_table_rhats <- function(res) {
    rhats <- rhat(res$fit, "foi")

    if (any(is.nan(rhats))) {
        rhats[which(is.nan(rhats))] <- 0
    }

    res_rhats <- data.frame(year = res$RealYexpo, rhat = rhats)
    # This is because I'm not estimating these foi values
    res_rhats$rhat[res_rhats$rhat == 0] <- NA

    return(res_rhats)
}

get_prev_expanded <- function(foi, dat) {
    ndat <- data.frame(age = 1:90)
    dim_foi <- dim(foi)[2]
    if (dim_foi < 90) {
        oldest_year <- 90 - dim_foi + 1
        foin <- matrix(NA, nrow = dim(foi)[1], 90)
        foin[, oldest_year:90] <- foi
        foin[, 1:(oldest_year - 1)] <- rowMeans(foi[, 1:5])
    } else {
        foin <- foi
    }

    foi_expanded <- foin


    age_class <- 1:NCOL(foi_expanded)
    ly <- NCOL(foi_expanded)
    exposure <- matrix(0, nrow = length(age_class), ncol = ly)
    for (k in 1:length(age_class)) exposure[k, (ly - age_class[k] + 1):ly] <- 1
    exposure_expanded <- exposure


    iterf <- NROW(foi_expanded)
    age_max <- NROW(exposure_expanded)
    PrevPn <- matrix(NA, nrow = iterf, ncol = age_max)
    for (i in 1:iterf) {
        PrevPn[i, ] <- 1 - exp(-exposure_expanded %*% foi_expanded[i, ])
    }

    lower <- apply(PrevPn, 2, function(x) quantile(x, 0.1))
    upper <- apply(PrevPn, 2, function(x) quantile(x, 0.9))
    medianv <- apply(PrevPn, 2, function(x) quantile(x, 0.5))


    predicted_prev <- data.frame(
        age = 1:90,
        predicted_prev = medianv,
        predicted_prev_lower = lower,
        predicted_prev_upper = upper
    )


    observed_prev <- dat %>%
        select(age_mean_f, frac, upperCI95, lowerCI95, N_seropos, N_tot) %>%
        rename(age = age_mean_f, sample_by_age = N_tot, positives = N_seropos)

    prev_expanded <- merge(predicted_prev, observed_prev, by = "age", all.x = TRUE) %>%
        mutate(survey = "")

    # I added this here for those cases when binned is prefered for plotting
    if (dat$Age_max[1] - dat$Age_min[1] < 3) {
        dat$cut_ages <- cut(as.numeric(dat$age_mean_f), seq(1, 90, by = 5),
            include.lowest = TRUE
        )
        xx <- dat %>%
            group_by(cut_ages) %>%
            summarise(
                bin_size = sum(total),
                bin_pos = sum(counts)
            )
        labs <- read.table(
            text = gsub("[^.0-9]", " ", levels(xx$cut_ages)),
            col.names = c("lower", "upper")
        ) %>%
            mutate(lev = levels(xx$cut_ages), midAge = round((lower + upper) / 2)) %>%
            select(midAge, lev)
        xx$midAge <- labs$midAge[labs$lev %in% xx$cut_ages]
        conf <- data.frame(Hmisc::binconf(xx$bin_pos, xx$bin_size, method = "exact"))
        xx <- cbind(xx, conf) %>% rename(
            age = midAge, p_obs_bin = PointEst,
            p_obs_bin_l = Lower, p_obs_bin_u = Upper
        )

        prev_final <- merge(prev_expanded, xx, by = "age", all.x = TRUE)
    } else {
        prev_final <- prev_expanded %>% mutate(
            cut_ages = "original",
            bin_size = sample_by_age,
            bin_pos = positives,
            p_obs_bin = frac,
            p_obs_bin_l = lowerCI95,
            p_obs_bin_u = upperCI95
        )
    }
    return(prev_final)
}



plot_info_table <- function(info, size_text) {
    dato <- data.frame(
        y = NROW(info):1,
        text = paste0(rownames(info), ": ", info[, 1])
    )
    p <- ggplot(dato, aes(x = 1, y = y)) +
        scale_y_continuous(limits = c(0, NROW(info) + 1), breaks = NULL) +
        # scale_x_continuous(breaks=NULL) +
        theme_void() +
        geom_text(aes(label = text), size = size_text / 2.5, fontface = "bold")

    return(p)
}






fCombinedPlots <- function(res, dat, lambda_sim = NA, max_lambda = NA, size_text = 25) {
    if (is.character(res$fit) == FALSE) {
        if (class(res$fit@sim$samples) != "NULL") {
            summary_mod <- extract_summary_mod(res, dat)
            foi <- rstan::extract(res$fit, "foi", inc_warmup = FALSE)[[1]]


            prev_expanded <- get_prev_expanded(foi, dat)
            plot_prev <-
                ggplot(prev_expanded) +
                geom_ribbon(aes(
                    x = age, ymin = predicted_prev_lower,
                    ymax = predicted_prev_upper
                ), fill = "#c994c7") +
                geom_line(aes(x = age, y = predicted_prev), colour = "#7a0177") +
                geom_errorbar(aes(age,
                    ymin = p_obs_bin_l,
                    ymax = p_obs_bin_u
                ), width = 0.1) +
                geom_point(aes(age, p_obs_bin, size = bin_size),
                    fill = "#7a0177", colour = "black", shape = 21
                ) +
                theme_bw(size_text) +
                coord_cartesian(xlim = c(0, 80), ylim = c(0, 1)) +
                theme(legend.position = "none") +
                ylab("Sero-positivity") +
                xlab("Age")



            #-------- This bit is to get the actual length of the foi data
            foi_dat <- res$foi_cent_est
            if (!is.na(lambda_sim)) {
                lambda_mod_length <- NROW(foi_dat)
                lambda_sim_length <- length(lambda_sim)

                if (lambda_mod_length < lambda_sim_length) {
                    remove_x_values <- lambda_sim_length - lambda_mod_length
                    lambda_sim <- lambda_sim[-c(1:remove_x_values)]
                }

                foi_dat$simulated <- lambda_sim
            }

            #--------
            foi_dat$medianv[1] <- NA
            foi_dat$lower[1] <- NA
            foi_dat$upper[1] <- NA


            plot_foi <-
                ggplot(foi_dat) +
                geom_ribbon(aes(x = year, ymin = lower, ymax = upper),
                    fill = "#41b6c4", alpha = 0.5
                ) +
                geom_line(aes(x = year, y = medianv),
                    colour = "#253494", size = size_text / 8
                ) +
                theme_bw(size_text) +
                coord_cartesian(ylim = c(0, max_lambda)) +
                # theme(panel.grid.major = element_blank(),
                #       panel.grid.minor = element_blank(),
                #       panel.background = element_blank(),
                #       axis.line = element_line(colour = "black")) +
                ylab("Force-of-Infection") +
                xlab("Year")


            if (!is.na(lambda_sim)) {
                lambda_plot <= lambda_plot + geom_line(aes(x = year, y = simulated),
                    colour = "red", size = 1.5
                )
            }

            rhats <- get_table_rhats(res)

            plot_rhats <- ggplot(rhats, aes(year, rhat)) +
                geom_line(colour = "purple") +
                geom_point() +
                coord_cartesian(ylim = c(0.7, 2)) +
                geom_hline(yintercept = 1.1, colour = "blue", size = size_text / 12) +
                theme_bw(size_text) +
                ylab("Convergence (R^)")


            # browser()
            plot_data <- plot_info_table(t(summary_mod), size_text = size_text)


            plot_arrange <- grid.arrange(plot_data,
                plot_prev,
                plot_foi,
                plot_rhats,
                nrow = 4,
                heights = c(1.5, 1, 1, 1)
            )
            # dev.off()
            res_plots <- list(
                plots = list(
                    plot_data = plot_data,
                    plot_prev = plot_prev,
                    plot_foi = plot_foi,
                    plot_rhats = plot_rhats
                ),
                summary_mod = summary_mod,
                rhats = rhats,
                prev_expanded = prev_expanded
            )
        }
    } else {
        print("model did not run")
        res_plots <- NA
    }
    return(res_plots)
}



#########################
# helpers for antibody  #
# prevalence extraction #
#########################

# this is to suppress and capture warnings and errors 
# more gracefully 
get_prev_expanded_safe = safely(get_prev_expanded)


# util function to create survey region codes 
# as unique IDs
get_survey_codes <- function(all_regions){
    fit_res_errs <- map(all_regions, function(study_id_i){
        # subset data
        dat <- filter(df_all_regions, study_id == study_id_i) %>% 
            arrange(., age_mean_f) %>% 
            mutate(., birth_year = Year - age_mean_f)
        # generate code for saved files
        survey_code = paste(dat$Country[1], dat$Region[1], dat$Year[1], sep='_')
        return(data.frame(study_id=study_id_i, file_id=survey_code))
    }) %>% bind_rows
}


# helper iterates over code IDs and extracts necessary data
get_prevalence <- function(all_regions, df_all_regions, df_file_id_path) {
    df_prevs <- map(
        all_regions,
        function(study_id_i) {
            cat(paste("Extracting from: ______ ______ ", study_id_i, "\n", sep = ""))
            # subset data
            dat <- filter(df_all_regions, study_id == study_id_i) %>%
                arrange(., age_mean_f) %>%
                mutate(., birth_year = Year - age_mean_f)
            # get file paths and read in data
            filepath <- df_file_id_path$filepath[df_file_id_path$file_id == study_id_i]
            read_res <- get(load(filepath)) #
            # extract foi
            foi <- rstan::extract(res$fit, "foi", inc_warmup = FALSE)[[1]]
            # extract prevalence
            prev_expanded <- get_prev_expanded_safe(foi, dat)
            df_prev_expanded <- prev_expanded$result
            # avoid crashing by defaulting return to NULL
            if (!is.null(df_prev_expanded)) {
                df_prev_expanded$survey_code <- study_id_i
            } else {
                df_prev_expanded <- NULL
            }
            return(df_prev_expanded)
        },
        .progress = T
    ) %>%
        compact() %>% # drop NULL items from list
        bind_rows() # and bind all into one long df
}




####################################
#### HELPERS FOR ERROR HANDLING ####
####################################

# capture fit results, errors, and warnings 
# using the sampler func wrapped around
# purrr::quietly and purrr::safely
safe_quiet_fit <- safely(.f = quietly(.f = sampling))
# capture loo cross-validation messages
safe_quiet_loo <- safely(.f = quietly(.f = loo))
# capture plotting messages
safe_quiet_fCombinedPlots <- safely(.f = quietly(.f = fCombinedPlots))


# since safely and quietly are nested, this cleans up the output
# albeit in a cumbersome way
tidy_up <- function(.inpt) {
    toreturn <- list(
        result = .inpt$result$result,
        output = .inpt$result$output,
        warnings = ifelse(length(.inpt$result$warnings) == 0, "", .inpt$result$warnings),
        messages = ifelse(length(.inpt$result$messages) == 0, "", .inpt$result$messages),
        error = .inpt$error$message
    ) %>% map(function(.npt) {
        # null for empty strings
        if ((typeof(.npt) == "character") && (nchar(.npt) == 0)) {
            .npt <- NULL
        }
        return(.npt)
    })
    return(toreturn)
}



###############################################
#### LOG MAKE/SAVE (ERROR/MESSAGE/WARNING) ####
###############################################



# creates a unique id based on date-time
# returns string "YYMMDD_HHMMSS"
get_dt_id <- function() {
    id <- as.character(Sys.time()) %>%
        substring(3, 19) %>%
        str_replace_all(":", "") %>%
        str_replace_all("-", "") %>%
        str_replace_all(" ", "_")
    return(id)
}


# writes all errors, messages, and warnings 
# can also save an additional CSV which is convenient
#### but needs to be specified here 
write_log_json <- function(
    fit_res_errs, # output of silenced fitter (list with results, output, error, messages, etc)
    vars_to_sub = NULL, # which variables from fit_res_errs to subset
    dest_dir = NULL, # where to save log, make sure to end with "/"
    save_csv = FALSE) {
    # defaults
    if (is.null(vars_to_sub)) {
        vars_to_sub <- c("study_id", "error", "warnings", "messages", "loo_warnings", "plot_warnings", "time_mins")
    }
    if (is.null(dest_dir)) {
        dest_dir <- paste0(getwd(), "/")
    }
    # loop over msgs for all regions
    df_log <- map(fit_res_errs, function(.item) {
        # loop over all vars_to_sub
        map(vars_to_sub, function(.field) {
            if (!is.null(.item[[.field]])) {
                # convert to df and name accordingly
                out <- data.frame(.x = paste0(.item[[.field]]))
                colnames(out) <- paste(.field)
                return(out)
            } else {
                NULL
            }
        }) %>% # drop nulls and bind cols together
            compact() %>% bind_cols() # data frame (ndim = 1 row and n cols)
    }) %>% bind_rows() # or ldply
    # save to dest_dir
    id <- get_dt_id()
    write(toJSON(df_log), file.path(dest_dir, paste0("log_", id, ".json")))
    if (save_csv == TRUE) {
        write_csv(df_log, file.path(dest_dir, paste0("log_", id, ".csv")))
    }
}



# dump all errors, warnings, and messages into a txt file
write_log_txt <- function(
    fit_res_errs, # output of silenced fitter (list with results, output, error, messages, etc)
    vars_to_sub = NULL, # which variables from fit_res_errs to subset
    dest_dir = NULL # where to save log, make sure to end with "/"
    ) {
    # defaults
    if (is.null(vars_to_sub)) {
        vars_to_sub <- c("region", "error", "warnings", "messages", "loo_warnings")
    }
    if (is.null(dest_dir)) {
        dest_dir <- getwd()
    }
    # create file
    id <- get_dt_id()
    fname <- file.path(dest_dir, paste0("log_", id, ".txt"))
    file.create(fname)
    # writing only warnings
    # map(fit_res_errs, c('warnings')) %>% unlist %>% cat(., file=fcon, append=T)
    map(fit_res_errs, function(.item) {
        map(vars_to_sub, function(.field) {
            if (!is.null(.item[[.field]])) {
                paste0(
                    toupper(.field), "\n", .item[[.field]], "\n\n"
                )
            } else {
                NULL
            }
        }) %>%
            compact() %>%
            unlist() # drop NULLs and convert to vector
    }) %>% # write file
        unlist() %>% cat(., file = fname, append = T, sep = "\n")
}

