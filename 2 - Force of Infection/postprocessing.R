library("purrr")
library("magrittr")
library("stringr")
library("tidyverse")
library("jsonlite")

# creating dirs and listing files 
data_dir <- file.path('res_kang/postprocessed')
if (!dir.exists(data_dir)) dir.create(data_dir)
res_dir = 'res_kang'
all_RData = list.files(res_dir, pattern = 'RData') # all files

# view logs and any error messages 
# e.g. by opening json file in vscode 
#### reformat using ctrl+A > ctrl+K > F
# (!!!) set log path 
log_path <- "res_kang/log_240416_153752.json"
json_logs <- read_json(log_path)
df_logs <- map(json_logs, as.data.frame) %>% bind_rows
# write.csv(df_logs, file='log_csv.csv', row.names=F)

# PROCESSING LOGS TO UNDERSTAND 
# ERRORS AND WARNINGS 
### rm links etc
df_logs$warnings = df_logs$warnings %>% 
    str_remove_all(.,'There were ') %>%
    str_remove_all(., 
        '. See\nhttps://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup\nto find out why this is a problem and how to eliminate them.'
) %>% str_remove_all(., 
    '. See\nhttps://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded'
    ) %>%
    str_remove_all(., '. See\nhttps://mc-stan.org/misc/warnings.html#bfmi-low')
# rm links etc 
df_logs$loo_warnings = df_logs$loo_warnings %>% 
    str_remove_all(., ". See help\\('pareto-k-diagnostic'\\) for details.\n") %>%
    str_replace_all(., 'slightly', 'SLIGHTLY') %>%
    str_replace_all(., 'too', 'TOO')
# sort by number of divergent transitions or chains with problems 
df_logs$warning_n = df_logs$warnings %>% str_split(., ' ') %>% sapply(., function(.x) .x[1]) %>% as.integer
df_logs = select(df_logs, -messages) %>% arrange(warning_n) %>% select(-warning_n)

# CHECKING FAILED PLOTS 
# indicate whether plotting was successful and make IDs 
no_plot = list.files('res_kang/fail') %>% str_remove_all('.png')
all_ids = list.files('res_kang', 'RData') %>% str_remove_all('.RData')


## NOT RUN WITH KANG ET AL 
# get country and year 
regions_fromid = map(all_ids, function(.entr) str_split_1(.entr, '_') %>% 
    .[-c(1, length(.))] %>% paste0(., collapse='_')) %>% unlist
countries_fromid = map(all_ids, function(.entr) str_split_1(.entr, '_') %>% .[1]) %>% unlist
years_fromid = map(all_ids, function(.entr) str_split_1(.entr, '_') %>% .[length(.)]) %>% unlist
fromid = data.frame(study_id = all_ids, region = regions_fromid, country = countries_fromid, year = years_fromid)
# bring all data frames together 
df_logs_full = left_join(df_logs, fromid, by='study_id') %>% select(study_id, region, country, year, warnings, loo_warnings)
df_logs_full$plotted = !(df_logs_full$region_id %in% no_plot)


# count types of warnings
df_logs <- df_logs %>% mutate(warning_type = case_when(
    str_detect(warnings, 'divergent transitions') ~ 'divergent transitions',
    str_detect(warnings, 'exceeded the maximum treedepth') ~ 'exceeded the maximum treedepth',
    str_detect(warnings, 'Bayesian Fraction of Missing Information') ~ 'Bayesian Fraction of Missing Information'
))

# save 
## save df_logs if study_id is present
write.csv(df_logs_full, file=file.path(data_dir, 'log_processed.csv'), row.names=F)



# check relative number of warnings 
df_logs$loo_warnings %>% table
df_logs$warnings %>% str_detect('divergent transitions') %>% sum(na.rm=T)
df_logs$warnings %>% str_detect('exceeded the maximum treedepth') %>% sum(na.rm=T)
df_logs$warnings %>% str_detect('Bayesian Fraction of Missing Information') %>% sum(na.rm=T)

# the same as above but in one line 
df_logs$warning_type %>% table



##################
# summarise fits # 
##################

# read in data to get region country etc 
data_path <- "data/all_data_kang.RData"
df_all_regions <- get(load(data_path))


# get country, region, and survey year from file names 
# bc i was not very bright and did not save that information in `res`
info_df = data.frame( # collapse to df and add surv_name
        study_id = str_remove(all_RData, '.RData'),  # for compatibility with logs above 
        surv_name = all_RData # with extension
        )

info_df <- df_all_regions %>% select(study_id, Country,Region, Year) %>% distinct %>% 
    left_join(info_df, by=join_by(study_id))      

colnames(info_df) %>% cat(sep=', ')

# pmap is much like apply(df, 1, fun(x)) i.e. apply by row 
# returns list of dfs
# takes a few seconds to run 
all_param_summary = purrr::pmap(info_df, function(study_id, Country, Region, Year, surv_name){
    example = get(load(file.path(res_dir, surv_name))) # load that survey res file
    # create a df containing all fois and identical country, region, study year cols
    example$foi_cent_est %>% mutate(
        Country=Country, Region=Region, Study_year=Year, study_id=study_id
        )
}, .progress = T)

# collapse list into one large df 
all_param_summary_df = all_param_summary %>% list_rbind()
# save(all_param_summary_df, file = 'all_param_summary.RData')


# summary of foi years available 
# repeats above for readability 
foi_years_summary = purrr::pmap(info_df, function(study_id, Country, Region, Year, surv_name){
    example = get(load(file.path(res_dir, surv_name))) # load that survey res file
    # create a df containing foi range as year_min-year_max
    example$foi_cent_est$year %>% range %>% paste(., collapse=' - ') %>%
        data.frame(
            Country=Country, Region=Region, Study_year=Year, study_id=study_id, foi_yrs=.
            )
}, .progress = T) %>%
    list_rbind()
# save(foi_years_summary, file = 'foi_years_summary.RData')


data_dir <- file.path('res_kang/postprocessed')
write.csv(all_param_summary_df, 
          file = file.path(data_dir, 'all_param_summary.csv'), row.names = F)
write.csv(foi_years_summary, 
          file = file.path(data_dir, 'foi_years_summary.csv'), row.names = F)




##########################
#### extract foi1000s ####
#### from posterior   ####
##########################

# country_region_year.RData 
file_ids = list.files('res_kang', 'RData') %>% str_remove_all('.RData')
filepaths = list.files('res_kang', 'RData', full.names = T)

# arrange by file size to gauge what the largest files to read in are
# (not essential) 
df_fpaths = data.frame(
    filepath = filepaths,
    size_MB = map(filepaths, function(.file) file.info(.file)$size / 1048576) %>% unlist
) %>% arrange(size_MB)

# iterate over file paths and IDs 
all_1000s = map2(filepaths, file_ids, function(.file, .id){
    # read in data file 
    read_res = get(load(.file))
    # extract the 1000 draws from the posterior 
    read_res$foi_post_1000s %>% as.data.frame %>% reshape2::melt(
        # reorganise into long format by year  
        value.name = 'foi_post_1000s', variable.name = 'year'
    ) %>% mutate(region_id = .id) 
    # since map returns lists, bind into one df and select relevant columns
}, .progress = T) %>% bind_rows %>% select(region_id, year, foi_post_1000s) 


# check and write file 
dim(all_1000s) 
write.csv(all_1000s, file=file.path(data_dir, 'foi_post_1000s.csv'), row.names=F)

# print available year ranges 
walk(file_ids, function(.id){
    dat = filter(all_1000s, region_id == .id)
    yrs = unique(dat$year) %>% droplevels %>% as.character %>% as.integer %>% sort
    cat(paste(paste0(range(yrs), collapse='-'), '    ', .id, '\n', sep='', collapse=''))
})




#####################
# extract antibody  # 
# prevalence levels # 
#####################

# file names are IDs based on survey region and year 
# read these in 
df_file_id_path <- data.frame(
    file_id = list.files("res_kang", "RData") %>% str_remove_all(".RData"),
    filepath = list.files("res_kang", "RData", full.names = T)
)


# this large data file is created using preprocessing.R which accesses `open_data`
## the file below is the example which is not geolocated 
# 'burkina_faso_gabon/data/df_burkina_faso_gabon_2015.RData'
data_path <- "data/all_data_kang.RData"
df_all_regions <- get(load(data_path))

# regions to be iterated over
# have to be unique 
all_regions = df_all_regions$study_id %>% unique

# iterates over all data and code IDs and extracts necessary data 
df_prevalence_ext <- get_prevalence(all_regions, df_all_regions, df_file_id_path)

# check and save file 
dim(df_prevalence_ext) 
write.csv(df_prevalence_ext, file=file.path(data_dir, 'antibody_prevalence.csv'), row.names=F)


