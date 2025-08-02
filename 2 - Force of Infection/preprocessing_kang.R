library("purrr")
library("tidyverse")
library("readxl")

# load and clean data
data_dir = 'data/KangSeroprevalenceDataFiles' #serological_data  ## change if data is elsewhere
# only interested in files that start with Sero
all_filenames = list.files(data_dir, 'Sero')
# read in all data into large list 
# all columns = characters 
all_data <- purrr::map(
    all_filenames,
    function(x) {
        read_excel(
            file.path(data_dir, x), 
            range = cell_cols(c("A","B","C","D","E","F","G")),
            col_types = c(
                "text","text","text","text","text","text","text"
            )
        )
    }
)

all_df = list_rbind(all_data)  # merge all data and convert to a data frame 

all_df %>% select(Country, Region, Year) %>% distinct %>% print(n=30)

# a study in Iran has one month olds 
all_df %>% filter(Age_min=="1 month")
all_df %>% filter(Country=="Iran")
# these are assumed to be in the 1 y.o. group
all_df$Age_min <- ifelse(all_df$Age_min=="1 month", 1, all_df$Age_min)

# some years are given as a range e.g. 2013-2014
# the first year is taken in this case
all_df <- all_df %>%
    mutate(
        Year = ifelse(nchar(Year) > 4, substr(Year, 1, 4), Year)
    ) %>%
    mutate(
        Year = as.numeric(Year),
        Age_min = as.numeric(Age_min),
        Age_max = as.numeric(Age_max),
        N_seropos = as.numeric(N_seropos),
        N_tot = as.numeric(N_tot)
    )

# round down age group mean 
all_df$age_mean_f = apply(all_df[,4:5], 1, mean) %>% floor 

# get fraction seropositives
all_df$frac = all_df$N_seropos / all_df$N_tot

# some fractions yield NAs; drop these 
all_df <- all_df %>% drop_na(frac) 


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
save(all_df, file='data/all_data_kang.RData')