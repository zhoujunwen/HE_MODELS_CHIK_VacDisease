library('ggplot2')
library('purrr')
library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library('tidyverse')
library('reshape2')
library('bayesplot')


fileids = list.files('res', 'RData') %>% str_remove_all('.RData')
filepaths = list.files('res', 'RData', full.names = T)

fileids = list.files('res_debug', 'RData') %>% str_remove_all('.RData')
filepaths = list.files('res_debug', 'RData', full.names = T)

fpath = filepaths[str_detect(filepaths, 'Kumbo')]
fid = fileids[str_detect(fileids, 'Kumbo')]

read_res = get(load(fpath))
kumbo_fit = read_res$fit

colnames(get_sampler_params(kumbo_fit)[[1]])

pairs(chain=1, x=kumbo_fit, pars=c("P", "lp__"))
pairs(chain=1, x=kumbo_fit, pars=c("sigma", "lp__"))

# does't work 
# ggsave(filename = 'res_debug/Cameroon_Kumbo_2007_pairs.png', plot=p_P)
