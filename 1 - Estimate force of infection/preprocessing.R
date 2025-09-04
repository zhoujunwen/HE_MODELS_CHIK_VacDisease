library('purrr')
library('tidyverse')

# load and clean data
data_dir = '../../open_access' #serological_data  ## change if data is elsewhere
# only interested in files that start with Sero
all_filenames = list.files(data_dir, 'Sero')
# read in all data into large list 
all_data = purrr::map(all_filenames, 
                      function(x) read.csv(paste0(data_dir,'/',x),
                                           colClasses = c('character','character',
                                                          'character','integer','integer',
                                                          'integer','integer')))
all_df = list_rbind(all_data)  # merge all data and convert to a data frame 
all_df = all_df[!is.na(all_df$Country),]  # two rows containing NAs removed
all_df = all_df[!is.na(all_df$N_seropos),]  # multiple empty N_seropos and N_tot removed  
all_df = rename(all_df,comments = X)

# assume upper age bound = 100 when not given (x3 rows)
all_df$Age_max[is.na(all_df$Age_max)] = 100 

# convert years to ints; subset first year for e.g. "2015-2016"
# will throw an NA coercion warning 
all_df$Year = ifelse(!is.na(as.integer(all_df$Year)),
                     all_df$Year, substr(all_df$Year,1,4)) %>% as.integer
# round down age group mean 
all_df$age_mean_f = apply(all_df[,4:5], 1, mean) %>% floor 

# get fraction seropositives
all_df$frac = all_df$N_seropos / all_df$N_tot

# binomial confidence intervals 
CIs = purrr::map2(all_df$N_seropos, all_df$N_tot, 
                  function(.x,.y) data.frame(t(
                      prop.test(x=.x, n=.y, 
                                conf.level=.95, 
                                correct=FALSE)$conf.int[1:2]))) %>% list_rbind(.)
colnames(CIs) = c('upperCI95', 'lowerCI95')
all_df = cbind(all_df, CIs)    # add to original data frame for plotting 

# create ID column 
all_df$study_id <- ifelse(
    is.na(all_df$Region), 
    paste(all_df$Country, all_df$Year),
    paste(all_df$Country, all_df$Region, all_df$Year)
    ) %>% str_replace_all(" ", "_")

# place it at the beginning 
all_df <- all_df %>% relocate(study_id, .before = Country)

# save data into file which will be used by fitting.R
save(all_df, file='data/all_data_open_access.RData')
