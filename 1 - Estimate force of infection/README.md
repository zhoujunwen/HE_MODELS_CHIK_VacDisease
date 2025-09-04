# Fitting a time-varying force of infection to CHIKV serological data

The following repository contains code which fits a time-varying force of infection (FOI) to open access seroprevalence data on CHIKV infection. The chosen model is compiled in Stan. There are three models: slow and fast epidemic models from the `serofoi` package and a mixture model with bimodal priors for FOI. 

The code has been adapted from [https://github.com/zmcucunuba/chfoi-col/tree/main](https://github.com/zmcucunuba/chfoi-col/tree/main) and from the [**serofoi** package](https://epiverse-trace.github.io/serofoi/index.html).

> **From local prevalence surveys to national burden of Chagas disease, a modelling pipeline. Authors: Julia Ledien, Zulma M. Cucunubá, Gabriel Parra-Henao, Eliana Rodríguez-Monguí, Andrew P. Dobson, Susana B. Adamo, María-Gloria Basáñe, Pierre Nouvellet.**


## Requirements 

The code was run using the following packages and their respective versions:


|   **package**   | **version** |
|:-----------:|:-------:|
|    Hmisc    |  5.1.0  |
|  bayesplot  | 1.10.0  |
|   cowplot   |  1.1.1  |
|    dplyr    |  1.1.3  |
|   epitrix   |  0.4.0  |
|   ggplot2   |  3.4.3  |
|    grid     |  4.3.1  |
|  gridExtra  |   2.3   |
|   gsubfn    |   0.7   |
|  jsonlite   |  1.8.7  |
|     loo     |  2.6.0  |
|  magrittr   |  2.0.3  |
|  parallel   |  4.3.1  |
|   pracma    |  2.4.2  |
|    purrr    |  1.0.2  |
|  reshape2   |  1.4.4  |
|    rstan    | 2.26.22 |
|   stringr   |  1.5.0  |
|  tidyverse  |  2.0.0  |
| vscDebugger (optional) |  0.5.2  |



## Geting started 

The code can be run in several simple steps, also reflected in the repo structure below: 

(1) Run preprocessing via `preprocessing.R`. This step may require additional manipulation if own data set is used. The data format is given by `data/all_data_open_access.RData`.

(2) Run loop fitting model to all regions (or unique region IDs) via `run.R`. 

(3) Run postprocessing script to collect and summarise all results via `postprocessing.R`.


> [!TIP]
> Run code via the `run.R` script. Some parameters cannot be changed directly there (e.g. the format of the log files is specified in `fitting.R`).

> [!IMPORTANT]  
> Ensure that the appropriate intialisation function (`utils::fInit()`) is chosen for the model specified in `run.R`. Also provide paths in the recommended format in all functions and files.


## Repo structure

If rerunning code ensure to follow the repo structure provided below and create a `./res/` dir. 

```
│
├───preprocessing.R     // (1) preprocess data before fitting (CIs etc)
├───run.R               // (2) read in data, fit models, and save results; 
│                       // sources `utils.R` and `fitting.R`
├───postprocessing.R    // (3) aggregate all FOIs together, 
├───utils.R             // adapted from original authors 
├───fitting.R           // large for-loop with to fit and save logs 
│                       // (mainly funcs for fitting/visualisation)
├───burkina_faso_gabon  // fitting the same model to data from Burkina Faso and Gabon 
│   ├───data            // (https://academic.oup.com/jid/article/227/2/261/
│   │                   // 6609553?login=false#supplementary-data)
│   └───res_bfg
├───data                // open access seroprevalence data (from Salje et al.)
├───mod                 // Stan models (choose appropriate model initialisation too)
└───res                 // RData files containing results (to be created)
    └───figs            // figures of (most) results 
```
