This repository contains the replication code for Roubenoff, Feehan, and Mahmud,
_Evaluating primary and booster vaccination prioritization strategies for COVID-19 by age and high-contact employment status using data from contact surveys_.


# Replication code and data files
This respository contains two R scripts containing the functions
requried for running the model and one RMarkdown file
containing the analysis. As well, there are three empty folders that will be required
for the model to run.

* `SEIR.R`: Contains the model (`SEI2RDV_ode(...)`), vaccination distribution (`vaccination(...)`),
and a wrapper function `sim(...)` that runs all 6 distribution scenarios for 
a supplied set of parameters. The latter function is used throughout the analysis. It returns
a list of two objects; a summary table `out` and a list of trajectories for each simulation, `trajectories`.
* `utils.R`: Contains utilities for generating the contact matrix (`employment_contact_matrix(..)`)
and functions for loading data and processing results.
* `SEIR_figures.Rmd`: Contains the analysis and generates figures.
* `data/`: Directory containing the replication data. You must download the data from the Harvard Dataverse; see below for details.
* `figs/`: Output directory for the figures used in the paper
* `sims_data/`: Cached simulations

We have deposited our data in the Harvard Dataverse (LINK). 
There are four files present in the repository:

* `national_wave4.csv`: BICS wave 4 respondents (selected columns only)
* `national_alters_wave4.csv`: BICS wave 4 respondents reported contacts (selected columns only)
* `acs_pop.csv`: Total population counts for children, adults, and seniors derived form the 2019 ACS 5-year estimates. This file was generated using the load_pop() function in utils.R.
* `jhu_data_raw.csv`: COVID-19 incidence data, retrieved from the [Johns Hopkins Data Repository]("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports_us/"). Retrieval of this file is outlined in utils.R.


# Replication procedure
To replicate our findings,
you will need to follow the following instructions. 

1) Clone this repository or download this replication code and place in your home 
directory. If you are unable to place this in that location, you must modify line 7
in SEIR.R (`setwd(...)`) to contain the location of the repository.
2) Download the four data files from the Harvard Dataverse and place them in 
the directory folder `data/`. 
3) Run all chunks in `SEIR_figures.Rmd`. It is not recommended that you 'knit' 
this file, rather just run all chunks in the file. 
> NB: This file can take a while to run (about 2 hours on a 2020 MacBook Air M1). We include a flag `run_sim` on line 28. When set to true (default), the script will cache the simulations in `sims_data/`. When set to false, the script will read the cached simulations without re-running them. Thus, we recommend running once with `run_sim<-TRUE`, and if additional runs are required, using `run_sim<-FALSE`.

