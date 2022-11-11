# SEIR utils

# Ethan Roubenoff


library(tidyverse)
library(reshape2)
library(deSolve)
library(tidycensus)
library(socialmixr)
library(cowplot)


#' Load ACS
#' 
#' @description Load the ACS pop and group into 3 child/adult/senior. Uses 2019 
#' ACS 5-year estimates of total population
#' 
#' @param none
#' @return data.frame 3 row x 2 column

load_pop <- function() {
  
  if (file.exists("data/acs_pop.csv")) {
    pop  <- read_csv("data/acs_pop.csv")
    return(pop)
  }
  
  tbl <- load_variables(2019, "acs5", cache = TRUE) 
  tbl <- tbl %>% 
    filter(name %in% paste0("B01001_0", sprintf("%02i", 1:49)) )
  
  pop <- get_acs("us", variables = tbl$name, cache = TRUE)
  pop <- pop %>% 
    left_join(tbl, by =c("variable"= "name"))
  pop <- pop %>% separate(label, sep = "!!|:!!|:", into = c("Estimate", "Total", "sex", "age"),
                          extra = "drop", fill = "right") %>%
                 mutate(age = na_if(age, ""),
                        sex = na_if(sex, ""))
  pop <- pop %>% 
    pivot_wider(id_cols = age, names_from = "sex", values_from = "estimate") %>%
    drop_na(age) 
  pop <- pop %>% select(-2)
  pop <- pop %>% mutate(acs_tot = Male + Female) %>% select(-Male, -Female)
  
  pop <- pop %>% mutate(ego_risk_category  = c(rep("child", 4), rep("adult", 13), rep("senior", 6))) %>%
    group_by(ego_risk_category) %>%
    summarize(acs_tot= sum(acs_tot))
  
  write_csv(pop, "data/acs_pop.csv")
  return(pop)
}

 





#' Load BICS respondents. Only works for wave 4.
#' 
#' @param wave=4 Wave to load
#' @param fix_unknowns=TRUE If TRUE, automatically redistribute
#' alters' unknown work status to match the ego proportion.

load_bics <- function(wave = 4, fix_unknowns = TRUE) {
  path <- "data/"
  
  if (wave == 3) {
    stop("Must be wave 4")
    national <- readRDS(file.path(path, "national_wave3.rds"))
    national_alters <- readRDS(file.path(path, "national_alters_wave3.rds"))
  } else if(wave == 4) {
    national <- read_csv(file.path(path, "national_wave4.csv"))
    national_alters <- read_csv(file.path(path, "national_alters_wave4.csv"))
  } else {
    stop("Invalid wave")
  }
  
  
  national <- national %>% mutate_at(vars(num_cc_work), ~if_else(. > 30, 30, .))
  
  agecat_remap <- function(x) addNA(cut(x, 
                                        breaks = c(0, 18, 65, 150), 
                                        labels = c("child", "adult", "senior"),
                                        include.lowest=TRUE,
                                        exclude=NULL,
                                        ordered=TRUE,
                                        right=FALSE))
  
  national <- national %>%
    mutate(agecat_ = agecat_remap(age))
  
  national <- national %>%
    mutate(ego_risk_category = as.character(agecat_)) %>%
    mutate(ego_risk_category = case_when(
      ego_risk_category == "adult" & !is.na(num_cc_work) & num_cc_work > 0 ~ "adult_working_inperson",
      ego_risk_category == "adult" ~ "adult_notworking_inperson",
      TRUE ~ ego_risk_category
    )) %>%
    mutate(ego_risk_category = addNA(factor(ego_risk_category, levels = c("child", "adult_notworking_inperson", 
                                                                          "adult_working_inperson", "senior"),
                                            ordered = TRUE)))

  
  ego_risk_cat_summary <- national %>%
    group_by(ego_risk_category) %>%
    summarize(n = n(), weighted_prop = sum(weight_pooled))
  
  national_alters <- national_alters %>%
    mutate(id = row_number())
  
  national_alters <- national_alters %>%
    mutate(
      alter_risk_category = case_when(
        alter_agecat == "[0,18)" ~ "child",
        alter_agecat == "[65,100]" ~"senior",
        (rel_work | rel_egoclient | rel_altercustomer | loc_store |loc_work ) ~ "adult_working_inperson",
        TRUE ~ "adult_unknown" 
      )) %>%
    drop_na(alter_risk_category)
  
  alters_risk_cat_summary <- national_alters %>%
    group_by(alter_risk_category) %>%
    summarize(n = n(), weighted_prop = sum(weight_pooled))
  

  ## This section randomly re-assigns high/low contact status to alters who 
  ## have unknown status such that the proportion of alters with each 
  ## status matches the egos 
  
  if (fix_unknowns) {
    proportion_egos_working_in_person <- 
      (ego_risk_cat_summary %>% 
         filter(ego_risk_category == "adult_working_inperson") %>% 
         pull(weighted_prop)) / 
      (ego_risk_cat_summary %>% 
         filter(ego_risk_category %in% c("adult_working_inperson", "adult_notworking_inperson")) %>% 
         pull(weighted_prop) %>%
         sum())  
    
    number_adult_alters <- 
      (alters_risk_cat_summary %>% 
         filter(alter_risk_category %in% c("adult_working_inperson", "adult_unknown")) %>% 
         pull(n) %>% sum())  
    
    n_alters_working_in_person <- 
      (alters_risk_cat_summary %>% 
        filter(alter_risk_category== "adult_working_inperson") %>% 
        pull(n)) 
    
    n_new <- round(proportion_egos_working_in_person*number_adult_alters - n_alters_working_in_person)
    
    unknown_adult_alters <-  national_alters %>%
      filter(alter_risk_category == "adult_unknown") %>%
      pull(id)
    
    new_alters_working_in_person <- sample(unknown_adult_alters, n_new)
    
    national_alters[new_alters_working_in_person, "alter_risk_category"] <- "adult_working_inperson"
    
    national_alters[
      national_alters$alter_risk_category == "adult_unknown", 
      "alter_risk_category" ] <- "adult_notworking_inperson"
    
  } else {
    national_alters[
      national_alters$alter_risk_category == "adult_unknown", 
      "alter_risk_category"
    ] <- "adult_notworking_inperson"
  }

  national_alters$alter_risk_category <- factor(national_alters$alter_risk_category, 
                                                levels = levels(national$ego_risk_category),
                                                ordered = TRUE)
  alters_risk_cat_summary <- national_alters %>%
    group_by(alter_risk_category) %>%
    summarize(n = n(), weighted_prop = sum(weight_pooled))
  

  
  
  national_alters <- left_join(national_alters,
                                     national %>% select(rid, ego_weight_pooled = weight_pooled, ego_risk_category),
                                     by = "rid")
  
  

  
  return(list(national = national,
              national_alters= national_alters,
              ego_risk_cat_summary = ego_risk_cat_summary,
              alters_risk_cat_summary= alters_risk_cat_summary))
}


#' Gets the corresponding POLYMOD matrix
#' 
#' @details Wrapper around socialmixr::contact_matrix. Code adapted from 
#' Ayesha Mahmud.
#' 
#' @param polymod_country Country to use contact survey data. Defaults to UK.
#' @param polymod_age_limits Cutoffs for age bins. 
#' @param school_contacts If FALSE, removes school contacts (primarily 0-18 age group)
#' @return (polymod_age_limits)^2 Matrix
getPolymodMatrix <- function(polymod_country = "United Kingdom",
                             polymod_age_limits = c(0, 18, 25, 35, 45, 65), 
                             school_contacts = FALSE){
  
  data(polymod)
  
  if(school_contacts){
    poly <- socialmixr::contact_matrix(polymod, 
                                       countries = polymod_country, 
                                       age.limits = polymod_age_limits, symmetric = TRUE)
    polymod_mat <- poly$matr
    
  } else {
    data_part <- polymod$participants
    data_cnt <- polymod$contacts %>% filter(cnt_school == 0)
    poly_no_school <- socialmixr::contact_matrix(survey(data_part, data_cnt), 
                                                 countries = polymod_country, 
                                                 age.limits = polymod_age_limits, symmetric = TRUE)
    polymod_mat <- poly_no_school$matrix
    
  }
  
  dimnames(polymod_mat)[[1]] <- dimnames(polymod_mat)[[2]]
  
  return(polymod_mat)
}




#' Calculates a contact matrix from Ego-Alter data
#' 
#' @details Adapted from Ayesha Mahmud
#' 
#' @param ego_df BICS respondents
#' @param alter_df BICS respondents' contacts
#' @param age_margins Total counts per age group (from ACS)
#' @param weightvar, alter_weightvar Name of weight variables in ego and alter dfs
#' @param agevar Name of age/risk category variable
#' @param wave0 ignored
#' 
#' @return Matrix of dimensions (length(age_margins))^2

calculate_contact_matrix <- function(ego_df, 
                                     alter_df, 
                                     age_margins,
                                     weightvar, 
                                     alter_weightvar=NULL, 
                                     agevar=NULL,
                                     wave0=FALSE) {
  
  
  # Rename vars
  wmap <- c('.weight'=weightvar)
  alter_wmap <- c('.ego_weight'=paste0("ego_",weightvar))
  
  ego_df <- ego_df %>%  rename(!!wmap)
  alter_df <- alter_df %>% rename(!!alter_wmap)
  
  if(is.null(alter_weightvar)) {
    alter_weightvar <- 'alter_weight' 
  }
  
  altervar_wmap <- c('.alter_weight'=alter_weightvar)
  alter_df <- alter_df %>% rename(!!altervar_wmap)
  
  if (! is.null(agevar)) {
    agemap <- c('.agecat'=paste0('ego_', agevar))
    alter_agemap <- c('.ego_agecat'=paste0('ego_', agevar),
                      '.alter_agecat'=paste0('alter_', agevar))
    
    ego_df <- ego_df %>%  rename(!!agemap)
    alter_df <- alter_df %>% rename(!!alter_agemap)
    age_margins <- age_margins %>% rename(!!agemap)
    
    if(! isTRUE(all.equal(levels(ego_df$.agecat),
                          levels(alter_df$.ego_agecat),
                          levels(alter_df$.alter_agecat),
                          levels(age_margins$.agecat)))) {
      
      stop('Levels for the age variable need to be the same in ego_df, alter_df, and targets\n.')
    }
  }
  
  
  ## get the denominator -- ie, the weighted number of interviews in each age group
  ## (need to do this separately b/c people who report 0 contacts in a given age group still
  ##  should be in the denominator when calculating avg number of connections to that age gp)
  
  mix_denom <- ego_df %>%
    select(.ego_age = .agecat, .weight) %>%     
    group_by(.ego_age) %>%                    
    summarize(num_interviews = n(),
              weighted_num_interviews = sum(.weight),
              .groups='drop_last') %>% 
    mutate(.ego_age = as.factor(.ego_age))
  
  unadj_contact_mat <- alter_df %>%     
    select(.ego_age=.ego_agecat, 
           .alter_age=.alter_agecat, 
           .alter_weight,          
           .ego_weight) %>%                                   
    group_by(.ego_age, .alter_age, .drop=FALSE ) %>%      
    # weighted_n_contacts is the weighted total number of reported contacts by respondents 
    #                     in ego_age to contacts in alter_age
    #      raw_n_contacts is the unweighted count of reported contacts by respondents
    #                     in ego_age to contacts in alter_age
    summarize(weighted_n_contacts = sum(.alter_weight*.ego_weight),
              raw_n_contacts = n(),
              .groups='drop_last') %>%
    # join on the denominator, ie, the number of interviews and weighted number of interviews
    # in ego_age
    left_join(mix_denom, by=c('.ego_age')) %>%       
    # join on the ego agegp size and alter agegp size in the population (from the ACS)
    left_join(age_margins %>%  
                select(.ego_age = .agecat, 
                       ego_acs_N = acs_tot),        
              by=c('.ego_age')) %>%                
    left_join(age_margins %>% 
                select(.alter_age = .agecat, 
                       alter_acs_N = acs_tot),    
              by=c('.alter_age')) %>%            
    # divide the populations by a million to keep the arithmetic from overflowing
    # (this will not affect the adjustment, as long as we divide ego_acs_N and alter_acs_N by
    #  the same thing)
    mutate(ego_acs_N = ego_acs_N / 1e6,
           alter_acs_N = alter_acs_N / 1e6)  %>%
    # unadj_avg_per_ego is the average number of contacts in age group alter_age
    # reported by respondents in ego_age, with the symmetry constraint not enforced
    mutate(unadj_avg_per_ego = weighted_n_contacts / weighted_num_interviews) %>%
    mutate(unadj_avg_per_ego = case_when(# if ego_acs_N is missing, this age group was not interviewed
      is.na(ego_acs_N) ~ NA_real_,
      # otherwise, this bootstrap resample may have happened not to resample
      # people in this age group - so the estimated number of contacts is 0
      is.na(weighted_num_interviews) ~ 0, 
      TRUE ~ unadj_avg_per_ego))
  
  
  # calculate the adjusted matrix, with symemtrization enforced
  # this symmetrization means that the matrix implies that
  #  n_e * avg_e,a = n_a * avg_a,e
  # where n_e is population total in ego age group
  #       n_a is population total in alter age group
  #       avg_e,a is avg number of connections to someone in a among people in e
  #       avg_a,e is avg number of connections to someone in e among people in a
  sym_contact_mat <- unadj_contact_mat %>%
    left_join( unadj_contact_mat %>% select(.alter_age = .ego_age, 
                                            .ego_age = .alter_age, 
                                            other_unadj_avg_per_ego = unadj_avg_per_ego,
                                            other_num_interviews = num_interviews, 
                                            other_weighted_num_interviews = weighted_num_interviews),     
               by=c('.ego_age', '.alter_age')) %>%                             
    mutate(sym_avg_per_ego = case_when(# if ego_acs_N is missing, this age group was not interviewed
      is.na(ego_acs_N) ~ NA_real_,
      # if alter_acs_N is missing, then the alter's age group was not interviewed 
      # (so no symmetrization is necessary)
      is.na(alter_acs_N) ~ unadj_avg_per_ego,
      # if alter's age group never reported contacts with ego age group (other_unadj_avg_per_ego will be NA)  
      # but ego's age group has reported contacts with alter, use that value
      is.na(other_unadj_avg_per_ego) & !is.na(alter_acs_N) ~ unadj_avg_per_ego,     ## added
      TRUE ~ (1 / (2*ego_acs_N)) * 
        ((unadj_avg_per_ego*ego_acs_N) + 
           (other_unadj_avg_per_ego*alter_acs_N)))) %>%
    mutate(fix_value = ifelse(is.na(other_unadj_avg_per_ego) & !is.na(alter_acs_N), 1, 0))    ## added
  
  sym_contact_mat <- ungroup(sym_contact_mat)
  
  return(sym_contact_mat)
}








#' Calculate employment contact matrix
#' 
#' @details explained in the supplementary matrieals
#' Calculated in the following steps:
#' 
#' Begin with the raw contact data. Derive three matrices:
#'   - Matrix A: BICS-derived work-status disaggregated (missing child contacts)
#'   - Matrix B: BICS-derived non-disaggregated (also missing child contacts)
#'   - Matrix C: POLYMOD-derived non-disaggregated
#' 
#' Process:
#'   - Scale matrix C to match dominant eigenvalue for matrix B
#'   - Use child-child entries for scaled C to fill in empty for mA
#' 
#' @param pop_data from load_pop() above
#' @param wave Survey wave passed to load_bics; must be 4
#' @param fix_unknowns Passed to load_bics
#' @return 4x4 Matrix 
#' 
#' 
employment_contact_matrix <- function(pop_data, wave = 4, fix_unknowns=TRUE, school_contacts=FALSE) {
  bics_data <- load_bics(wave = wave, fix_unknowns = fix_unknowns)
  agecat_remap <- function(x) addNA(cut(x,
                                        breaks = c(0, 18, 65, 150),
                                        labels = c("child", "adult", "senior"),
                                        include.lowest=TRUE,
                                        exclude=NULL,
                                        ordered=TRUE,
                                        right=FALSE))
  
  
  #' # Matrix A: BICS-derived work-status disaggregated (missing child contacts)
  #' 
  #' Prep data
  ## --------------------------------------------------------------------------------------
  # Scale the pop data margins with the survey proportion for in-person work
  
  pop_data_employment <- pop_data
  
  # Split adults into working in person/not working in person based on weighted
  # proportions from the survey
  proportions <- bics_data$ego_risk_cat_summary %>% 
    pull(weighted_prop)
  names(proportions) <- bics_data$ego_risk_cat_summary$ego_risk_category
  proportions <- proportions[names(proportions) != "senior"] 
  proportions <- proportions/sum(proportions)
  adult_working_inperson_proportion <- proportions["adult_working_inperson"]
  adult_notworking_inperson_proportion <- proportions["adult_notworking_inperson"]
  pop_data_employment[4,] <- list("adult_working_inperson", 200484607 * adult_working_inperson_proportion) 
  pop_data_employment[5,] <- list("adult_notworking_inperson", 200484607 * adult_notworking_inperson_proportion) 
  pop_data_employment <- pop_data_employment[-1,]
  
  row_order <- c("child","adult_notworking_inperson","adult_working_inperson","senior")
  
  
  #' Calculate matrix A 
  ## ---------------------------------------------------------------------------
  mA <- calculate_contact_matrix(ego_df = bics_data$national,
                                 alter_df = bics_data$national_alters,
                                 age_margins = pop_data_employment,
                                 weightvar = "weight_pooled",
                                 alter_weightvar = "alter_weight",
                                 agevar = "risk_category")
  
  mA <- mA %>% drop_na(.ego_age, .alter_age)
  
  mA <- mA %>% select(.ego_age, .alter_age, sym_avg_per_ego) 
  mA <- acast(mA, .ego_age ~ .alter_age, value.var = "sym_avg_per_ego")
  mA <- mA[row_order,row_order]
  
  
  #' # Matrix B: BICS-derived non-disaggregated (also missing child contacts)
  ## ---------------------------------------------------------------------------
  mB <- calculate_contact_matrix(ego_df = bics_data$national %>% 
                                   mutate(ego_risk_category= agecat_remap(age)),
                                 alter_df = bics_data$national_alters %>% 
                                   mutate(alter_risk_category= agecat_remap(alter_age),
                                          ego_risk_category = agecat_remap(ego_age)),
                                 age_margins = pop_data,
                                 weightvar = "weight_pooled",
                                 alter_weightvar = "alter_weight",
                                 agevar = "risk_category")
  
  mB <- mB %>%
    drop_na(.ego_age, .alter_age)
  
  mB <- mB %>% select(.ego_age, .alter_age, sym_avg_per_ego) 
  mB <- acast(mB, .ego_age ~ .alter_age, value.var = "sym_avg_per_ego")
  mB <- mB[c("child", "adult", "senior"), c("child", "adult", "senior")]
  
  
  #' # Matrix C: POLYMOD-derived non-disaggregated
  ## ---------------------------------------------------------------------------
  mC <- getPolymodMatrix(polymod_country = "United Kingdom", 
                         polymod_age_limits = c(0, 18, 65),
                         school_contacts = school_contacts)
  
  colnames(mC) <- rownames(mC) <- c("child", "adult", "senior")
  
  
  #' Calculate ratio of eigenvalues
  #' Contact in BICS is about half what it is in POLYMOD
  ## ---------------------------------------------------------------------------

  ratio <- eigen(mB)$values[1] / eigen(mC)$values[1]
  ratio
  
  #' Fill mA with child vallues from mC
  ## ---------------------------------------------------------------------------
  mA[1,1] <- mC[1,1] * ratio
  
  return(mA)
}




#' Processes the trajectories into data frames

process_SEI2RDV_trajectory <- function(trajectory) {
  vec <- colnames(trajectory)
  vec <- vec %>%
    str_replace("1", "_child") %>%
    str_replace("2", "_adult_notworking_inperson") %>%
    str_replace("3", "_adult_working_inperson") %>%
    str_replace("4", "_senior")
  
  colnames(trajectory) <- vec
  
  trajectory <- as_tibble(as.data.frame(trajectory))
  
  return(trajectory)
}



# To get the latest JHU data
github_url <- paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports_us/", 
                     "01-01-2021",
                     ".csv")
if (!file.exists("data/jhu_data_raw.csv")) {
  jhu_data_raw <- read_csv(github_url,col_types = cols())
  write_csv(jhu_data_raw, "data/jhu_data_raw.csv")
} else {
  jhu_data_raw <- read_csv("data/jhu_data_raw.csv")
}

#' Create the starting population
#' 
#' @details Pulls COVID rates from the Johns Hopkins respository
#' 
#' @param pop from load_pop
#' @param rho 
create_SEI2RDV_starting_pop <- function(pop, 
                                        start_date = "01-01-2021", 
                                        rho = c(0.35, 0.4, 0.4, 0.75), 
                                        u = c(0.39, 0.83, 0.83, 0.74),
                                        mu = c(0, 0.001, 0.001, 0.01),
                                        adult_prop = c(adult_notworking_inperson = 0.66, 
                                                       adult_working_inperson = 0.34),
                                        latentPeriod = 3,
                                        infectiousPeriod = 5,
                                        adjust = FALSE) {
  
  jhu_data <- jhu_data_raw %>%
    filter(ISO3 == "USA") %>%
    filter(!Province_State %in% c("Diamond Princess", "Grand Princess")) %>%
    select(Province_State, Confirmed, Deaths, Active, Recovered)
  
  # For states that have "Active" and "Recovered" listed, take proportions of total confirmed:
  counts <- jhu_data %>%
    select(Confirmed, Deaths, Recovered, Active) %>%
    drop_na() %>%
    colSums() 
  proportions <- counts / counts[1]
  
  # Scale total confirmed by these proportions
  counts <- jhu_data %>% 
    select(Confirmed, Deaths, Recovered, Active)  %>%
    colSums()
  counts["Recovered"] <- counts["Confirmed"] * proportions["Recovered"]
  counts["Active"] <- counts["Confirmed"] * proportions["Active"]
    
  
  
  
  
  # Create init matrix
  categories <- c("child", "adult_notworking_inperson", 
                  "adult_working_inperson", "senior")
  compartments <- c("N", "S", "Sx", "E", "Ic", "Isc", "R", "D", "Va", "Vb", "Vbx",
                    "Vboost", "Vc", "Cc", "Csc") 
  init <- matrix(0, nrow = length(compartments), ncol = 4, 
                 dimnames = list(compartments,  
                                 categories))
  
  # Create pop
  n_adults <- pull(pop[pop$ego_risk_category == "adult", "acs_tot"])
  pop <- bind_rows(pop, tibble(ego_risk_category = names(adult_prop), 
                               acs_tot = adult_prop * n_adults))
  pop <- pop %>% filter(ego_risk_category!="adult")
  pop_list <-pop$acs_tot
  names(pop_list) <- pop$ego_risk_category
  
  init["N", names(pop_list)] <- pop_list
  init["Ic", ] <- counts["Active"] * (u * init["N",])/sum(u * init["N",])
  init["D", ]   <- c(0,0,0,0)
  init["R", ]   <- c(0,0,0,0)
  init["R",]    <- c(0,0,0,0)
  init["Isc", ] <- c(0,0,0,0)
  init["E", ]   <- (init["Ic", ]) * (latentPeriod/infectiousPeriod)
  

  # Set Susceptibles to be everyone else
  init["S", ] <- init["N", ] - colSums(init[-c(1), ])
  

  init <- round(init)
  return(init)
}




#### Methods for calculating transmission probs #####

compute_R0 = function(u,          # category-dependent susceptibility to infection 
                      C,          # contact matrix
                      gamma,      # recovery rate
                      rho,        # category-dependent probability of having clinical infections
                      alpha,      # relative infectiousness of subclinical vs clinical 
                      n = 4){     # number of categories
  
  C <- as.matrix(C)
  # Davies NGM
  Du <- diag(u, n)
  Dy <- diag(1/gamma, n)
  Dx <- diag((1-rho)*alpha+rho, n)
  NGM <- Du %*% C %*% Dy %*% Dx 
  R0  <- abs(eigen(NGM)$values[1])
  return(R0)
}



scale_u_for_R0 = function(wanted_R0 = 2.5,  # The R_0 to optimize beta for
                          # category-dependent susceptibility to infection (taylor calls this u)
                          u =  c(0.39, 0.8, 0.8, 0.74), 
                          # contact matrix
                          C,      
                          # recovery rate
                          gamma = 1/6,      
                          # category-dependent probability of having clinical infections
                          rho = c(0.4, 0.75, 0.75, 0.9),  
                          # relative infectiousness of subclinical vs clinical
                          alpha = 0.5){      
  ret <- optimize(
    function(scale, 
             target, 
             beta, 
             C, 
             gamma, 
             rho, 
             alpha) {
      abs(target - compute_R0(u= u/scale, C = C, gamma, rho = rho, alpha = alpha))
  },
    interval = c(0.1, 40),
    target = wanted_R0,
    beta = beta,
    C = C,
    gamma = gamma,
    rho = rho,
    alpha = alpha
  )
  
  return(ret$minimum)
  
}






#### Calculate Stats ####

calculate_SEI2RDV_stats <- function(trajectories,          # List of simulations
                                    params,                # Parameter list (used for text output)
                                    printout = TRUE,       # Include text output?
                                    stats_start_date = 1   # Calculate incidence and deaths starting x days after simulation start
) {
  stats_df <- tibble(sim = rep(NA_character_, length(trajectories)),
                     total_deaths = rep(NA_real_, length(trajectories)),
                     total_sick = rep(NA_real_, length(trajectories)),
                     total_clinical = rep(NA_real_, length(trajectories)),
                     total_vaccinated= rep(NA_real_, length(trajectories)),
                     highest_sick = rep(NA_real_, length(trajectories))
  )
  
  trajectories <- map(trajectories, ~process_SEI2RDV_trajectory(.))
  
  stats_df$sim <- names(trajectories)
  
  # For all, start on day stats_start_date
  if (is.na(stats_start_date)) {
    stats_start_date <- 0
  }
  
  # Total Deaths
  stats_df[, "total_deaths"] <- map_dbl(trajectories, ~{
    as.data.frame(.) %>%
      slice(., {{stats_start_date}}, n()) %>%
      select(starts_with("D_")) %>%
      rowSums(.) %>%
      {.[2] - .[1]}
  })
  # Total number infected
  stats_df[, "total_sick"] <- map_dbl(trajectories, ~{
    as.data.frame(.) %>%
      slice(., {{stats_start_date}}, n()) %>%
      select(starts_with("Cc_"), starts_with("Csc_")) %>%
      rowSums(.)%>%
      {.[2] - .[1]}
  })
  # Total number clinical 
  stats_df[, "total_clinical"] <- map_dbl(trajectories, ~{
    as.data.frame(.) %>%
      slice(., {{stats_start_date}}, n()) %>%
      select(starts_with("Cc_")) %>%
      rowSums(.) %>%
      {.[2] - .[1]}
  })
  # Total number vaccinated 
  stats_df[, "total_vaccinated"] <- map_dbl(trajectories, ~{
    as.data.frame(.) %>%
      slice(., n()) %>%
      select(starts_with("Vc_")) %>%
      rowSums(.)
  })
  # Maximum sick
  stats_df[, "highest_sick"] <- map_dbl(trajectories, ~{
    # as.data.frame(.)%>%
    #   rowwise() %>%
    #   mutate(total_sick = sum(c_across(starts_with("Ic_")))) %>%
    #   pull(total_sick) %>%
    #   max(.)
    max(rowSums(select(as.data.frame(.), starts_with("Ic_"))))
  })
  
  
  
  
  if (printout) {
    n_deaths_novax <- pull(stats_df[stats_df$sim == "trajectory_none", "total_deaths"])
    n_deaths_seniorvax <- pull(stats_df[stats_df$sim == "trajectory_senior", "total_deaths"])
    n_deaths_hr <- pull(stats_df[stats_df$sim == "trajectory_hr", "total_deaths"])
    n_deaths_split <- pull(stats_df[stats_df$sim == "trajectory_split", "total_deaths"])
    n_deaths_tiered_sr <- pull(stats_df[stats_df$sim == "trajectory_tiered_sr", "total_deaths"])
    n_deaths_tiered_hc <- pull(stats_df[stats_df$sim == "trajectory_tiered_hc", "total_deaths"])
    deaths_list <- list(
      novax = n_deaths_novax, 
      seniorvax = n_deaths_seniorvax, 
      hr = n_deaths_hr, 
      split = n_deaths_split,
      tiered_sr = n_deaths_tiered_sr,
      tiered_hc = n_deaths_tiered_hc
      )
    least_deaths <-  names(deaths_list)[which.min(deaths_list)]
    
    n_clinical_novax <- pull(stats_df[stats_df$sim == "trajectory_none", "total_clinical"])
    n_clinical_seniorvax <- pull(stats_df[stats_df$sim == "trajectory_senior", "total_clinical"])
    n_clinical_hr <- pull(stats_df[stats_df$sim == "trajectory_hr", "total_clinical"])
    n_clinical_split <- pull(stats_df[stats_df$sim == "trajectory_split", "total_clinical"])
    n_clinical_tiered_sr <- pull(stats_df[stats_df$sim == "trajectory_tiered_sr", "total_clinical"])
    n_clinical_tiered_hc <- pull(stats_df[stats_df$sim == "trajectory_tiered_hc", "total_clinical"])
    
    clinical_list <- list(
      novax = n_clinical_novax, 
      seniorvax = n_clinical_seniorvax, 
      hr = n_clinical_hr, 
      split = n_clinical_split,
      tiered_sr = n_clinical_tiered_sr,
      tiered_hc = n_clinical_tiered_hc
      )
    least_clinical <-  names(clinical_list)[which.min(clinical_list)]
    
    
    n_highest_sick_novax <- pull(stats_df[stats_df$sim == "trajectory_none", "highest_sick"])
    n_highest_sick_seniorvax <- pull(stats_df[stats_df$sim == "trajectory_senior", "highest_sick"])
    n_highest_sick_hr <- pull(stats_df[stats_df$sim == "trajectory_hr", "highest_sick"])
    n_highest_sick_split <- pull(stats_df[stats_df$sim == "trajectory_split", "highest_sick"])
    n_highest_sick_tiered_sr <- pull(stats_df[stats_df$sim == "trajectory_tiered_sr", "highest_sick"])
    n_highest_sick_tiered_hc <- pull(stats_df[stats_df$sim == "trajectory_tiered_hc", "highest_sick"])
    
    highest_sick_list <- list(
      novax = n_highest_sick_novax, 
      seniorvax = n_highest_sick_seniorvax, 
      hr = n_highest_sick_hr, 
      split = n_highest_sick_split,
      tiered_sr = n_highest_sick_tiered_sr,
      tiered_hc = n_highest_sick_tiered_hc)
    least_highest_sick <-  names(highest_sick_list)[which.min(highest_sick_list)]
    
    message("------------------------------- \n",
            "------Simulations Summary------ \n",
            "------------------------------- \n\n",
            "Simulation Parameters:",
            "\n Theoretical R0: ", params$R0,
            "\n Derived beta (per-contact probability of transmission): ", params$beta, 
            "\n Latent Period: ", 1/params$sigma,
            "\n Infectious Period: ", 1/params$gamma,
            "\n alpha (Relative transmissibility of subclinical vs. clinical): ", params$alpha,
            "\n rho (proportion of each age group clinicial): ", paste(params$rho, collapse = ", " ),
            "\n mu (mortality probability in clinical): ", paste(params$mu, collapse = ", "),
            "\n v0 (initial day of vaccination): ", params$v0,
            "\n Number of daily vaccinations: ", params$daily_vax,
            "\n Simulation start date: ", params$start_date,
            "\n\n",
            
            
            
            
            "Relative to vaccinating no one,\n vaccinating seniors saves ",
            round(n_deaths_novax - n_deaths_seniorvax), 
            " lives,\n a ",  round((n_deaths_novax - n_deaths_seniorvax)/n_deaths_novax * 100), "% reduction. \n",
            
            "Under this scenario,\n the number of clinical infections is reduced by ",
            round(n_clinical_novax - n_clinical_seniorvax), 
            ",\n a ",  round((n_clinical_novax - n_clinical_seniorvax)/n_clinical_novax * 100), "% reduction. \n\n",
            
            
            "Relative to vaccinating no one,\n vaccinating high risk workers saves ",
            round(n_deaths_novax - n_deaths_hr), 
            " lives,\n a ",  round((n_deaths_novax - n_deaths_hr)/n_deaths_novax * 100), "% reduction. \n",
            
            "Under this scenario,\n the number of clinical infections is reduced by ",
            round(n_clinical_novax - n_clinical_hr), 
            ",\n a ",  round((n_clinical_novax - n_clinical_hr)/n_clinical_novax * 100), "% reduction. \n\n",
            
            
            "Relative to vaccinating no one,\n splitting vaccines between seniors and high risk workers saves ",
            round(n_deaths_novax - n_deaths_split), 
            " lives,\n a ",  round((n_deaths_novax - n_deaths_split)/n_deaths_novax * 100), "% reduction. \n",
            
            "Under this scenario,\n the number of clinical infections is reduced by ",
            round(n_clinical_novax - n_clinical_split), 
            ",\n a ",  round((n_clinical_novax - n_clinical_split)/n_clinical_novax * 100), "% reduction. \n \n \n \n",
            
            
            
            "Relative to vaccinating no one,\n a tiered system that first ",
            "prioritizes seniors \n and then HC workers saves ",
            round(n_deaths_novax - n_deaths_tiered_sr), 
            " lives,\n a ",  round((n_deaths_novax - n_deaths_tiered_sr)/n_deaths_novax * 100), "% reduction. \n",
            
            "Under this scenario,\n the number of clinical infections is reduced by ",
            round(n_clinical_novax - n_clinical_tiered_sr), 
            ",\n a ",  round((n_clinical_novax - n_clinical_tiered_sr)/n_clinical_novax * 100), "% reduction. \n \n \n \n",
            
            
            
            "Relative to vaccinating no one,\n a tiered system that first ",
            "prioritizes HC workers \n and then seniors saves ",
            round(n_deaths_novax - n_deaths_tiered_hc), 
            " lives,\n a ",  round((n_deaths_novax - n_deaths_tiered_hc)/n_deaths_novax * 100), "% reduction. \n",
            
            "Under this scenario,\n the number of clinical infections is reduced by ",
            round(n_clinical_novax - n_clinical_tiered_hc), 
            ",\n a ",  round((n_clinical_novax - n_clinical_tiered_hc)/n_clinical_novax * 100), "% reduction. \n \n \n \n")
            
    
            nameremap <- c(
              "novax" = "No one",
              "seniorvax" = "65+",
              "hr" = "High Contact Workers",
              "split" = "Split 65+ and Workers",
              "tiered_sr" = "Tiered (65+ then Workers)",
              "tiered_hc" = "Tiered (Workers then 65+)"
            )
            
            
            deaths_ct <- map_dbl(deaths_list, ~round((.-deaths_list[[least_deaths]])))
            deaths_pct <- map_dbl(deaths_list, ~round((.-deaths_list[[least_deaths]]) / . * 100, 2))
            names(deaths_pct) <- nameremap
            clinical_ct <- map_dbl(clinical_list, ~round((.-clinical_list[[least_clinical]])))
            clinical_pct <- map_dbl(clinical_list, ~round((.-clinical_list[[least_clinical]]) / . * 100, 2))
            names(clinical_pct) <- nameremap
    
            message("The most lives are saved by vaccinating ", nameremap[least_deaths],
                    paste0("\n a ", deaths_pct, "% reduction (",deaths_ct, ") over prioritizing ", names(deaths_pct)),
                    "\n \n \n \n"
            )
            
            message("The most clinical infections are spared by vaccinating ", nameremap[least_clinical],
                    paste0("\n a ", clinical_pct, "% reduction (", clinical_ct, ") over prioritizing ", names(deaths_pct)),
                    "\n \n \n \n"
                    
            )
            
            
            
            message(
            "------------------------------- \n",
            "-----------End Summary--------- \n",
            "------------------------------- \n"
            
    )
    
    
  }
  
  
  return(stats_df)
}





# Function to aggregate the CDC data on vaccine uptake
load_vu <- function() {
  
  # Load variables
  tbl <- load_variables(2019, "acs5", cache = TRUE) 
  tbl <- tbl %>% 
    filter(name %in% paste0("B01001_0", sprintf("%02i", 1:49)) )
  
  # Retrieve population and reformat
  pop <- get_acs("us", variables = tbl$name, cache = TRUE)
  pop <- pop %>% 
    left_join(tbl, by =c("variable"= "name"))
  pop <- pop %>% separate(label, sep = "!!|:!!|:", into = c("Estimate", "Total", "sex", "age"),
                          extra = "drop", fill = "right") %>%
    mutate(age = na_if(age, ""),
           sex = na_if(sex, ""))
  pop <- pop %>% 
    pivot_wider(id_cols = age, names_from = "sex", values_from = "estimate") %>%
    drop_na(age) 
  pop <- pop %>% select(-2)
  pop <- pop %>% mutate(acs_tot = Male + Female) %>% select(-Male, -Female)
  
  # Group by age and summarize
  pop <- pop %>% mutate(age_group  = c(
    "Under 5",
    rep("5-17", 3), 
    rep("18-29", 5),
    rep("30-39", 2),
    rep("40-49", 2),
    rep("50-64", 4),
    rep("65-74", 3),
    rep("75+", 3))
    ) %>%
    mutate(age_group = as_factor(age_group)) %>%
    group_by(age_group) %>%
    summarize(acs_tot= sum(acs_tot))
  
  
  
  # https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Demographics-in-the-United-St/km4m-vcsb
  # Week of Apr25-May1 2021, Fully vaccinated with 1 or 2 doses
  vu_apr2021 <- pop %>%
    filter(age_group != "Under 5",
           age_group != "5-17") %>%
    bind_cols(vu = c(20.8, 25.8, 36.0, 49.0, 76.8, 81.6))
  
  # Collapse into just2 age groups
  vu_apr2021 <- vu_apr2021 %>% mutate(
    age_group = case_when(
      age_group == "65-74" | age_group == "75+" ~ "Senior",
      TRUE ~ "Adult"
    )
  )
  
  # Take weighted average
  vu_apr2021 <- vu_apr2021 %>% 
    mutate(age_group = as_factor(age_group)) %>%
    group_by(age_group) %>%
    summarize(vu = sum(vu*acs_tot)/sum(acs_tot)) 
  
  vu_apr2021 <- bind_rows(list(age_group = "Child", vu = 0L), vu_apr2021)
  
  
  # Repeat with January 9-15 2022 to figure out boosters
  vu_jan2022 <- pop %>%
    filter(age_group != "Under 5") %>%
    bind_cols(vu = c(34, 67.5, 76.7, 80.3, 87.8, 95.3, 94.4),
              booster = c(0, 36.7, 43.4, 49.8, 61.2, 70.1, 76.4)) %>%
    mutate(
      age_group = case_when(
        age_group == "5-17" ~ "Child",
        age_group == "65-74" | age_group == "75+" ~ "Senior",
        TRUE ~ "Adult"
      )) %>%
    mutate(age_group = as_factor(age_group)) %>%
    group_by(age_group) %>%
    summarize(vu = sum(vu*acs_tot)/sum(acs_tot),
              booster = sum(booster*acs_tot)/sum(acs_tot))
  
  # Need to split these into 4-element lists
  
  vu_apr2021 <- vu_apr2021 %>% deframe() %>% .[c(1,2,2,3)]/100
  booster <- vu_jan2022 %>% select(age_group, booster) %>% deframe() %>% .[c(1,2,2,3)]/100
  vu_jan2022 <- vu_jan2022 %>% select(age_group, vu) %>% deframe() %>% .[c(1,2,2,3)]/100

  return(list(vu_apr2021 = vu_apr2021,
              vu_jan2022 = vu_jan2022,
              booster = booster
              ))
}









# Turn a matrix into a named list, where names are of the form rowname_colname
listmat <- function(m) {
  l <- as.vector(t(m))
  
  names(l) <- expand_grid(v1 = rownames(m) , v2 = colnames(m)) %>%
    mutate(v3 = paste0(v1, "_", v2)) %>%
    pull(v3)
  
  return(l)
}

# Turns a named list into a matrix, where names are of the form rowname_colname
unlistmat <- function(l) {
  
  l_names <- str_split(names(l), "_", 2)
  row_names <- unique(map_chr(l_names, ~.[1]))
  col_names <- unique(map_chr(l_names, ~.[2]))
  
  m <- matrix(l, nrow = length(row_names), ncol = length(col_names), 
         dimnames = list(row_names, col_names),
         byrow = TRUE)
  
  return(m)
}



# New function figures out the number of people we want to vaccinate each day

vr_func2 <- function(
  vr, # scalar number of vaccinations per day
  vmax, # 4-element vector, number of people in each group waiting to be vaccinated
  priority_queue = list(), # List of groups in order of priority. 
                           # Special case: if priority is NA then no one gets vaccinated
  split_proportion = c("adult_working_inperson" = 0.5, "senior" = 0.5) 
                         
) {
  
  # Check input type
  if (typeof(priority_queue) != "list") {stop("Priority queue MUST be a list!")}
  
  # Initialize counts
  priority_vaccinations <- vmax
  priority_vaccinations[] <- 0 # Create 0-vector with appropriate names
  vr_remaining <- vr - sum(priority_vaccinations)
  
  # Check if we are in null scenario
  if (length(priority_queue) >= 1) {
    if (any(is.null(priority_queue[[1]])) | 
        any(is.na(priority_queue[[1]]))) {
      return(priority_vaccinations)
      }
  }
  
  
  for (i in priority_queue) {
    # Handle split proportion
    if (length(i) > 1) {vr_remaining <- vr_remaining*split_proportion[i]}
    priority_vaccinations[i] <- floor(priority_vaccinations[i] + pmax(pmin(vmax[i], vr_remaining), 0))
    
    # How many doses are left over
    vr_remaining <- vr - sum(priority_vaccinations)
    
    # How many are still awaiting vaccination
    vmax <- vmax - priority_vaccinations
    
    # Check to see there aren't any priority people still waiting (seems to be issue with split)
    if (vr_remaining < 1) return(priority_vaccinations)
    
    # Edge case with split vax: need to ensure that if the first group is split, any
    # leftover vaccines go to the other prioritized group before they are distributed out
    if (any(vmax[i] > 0)) {
      # Drop the one that is 0
      i <- i[-c(which(vmax[i] == 0))]
      
      # Re-distribute
      if (length(i) > 1) {vr_remaining <- vr_remaining*split_proportion[i]}
      priority_vaccinations[i] <- floor(priority_vaccinations[i] + 
                                          pmax(pmin(vmax[i], vr_remaining), 0))
      vr_remaining <- vr - sum(priority_vaccinations)
      vmax <- vmax - priority_vaccinations
    }
  }
  
  
  # Split remaining vaccines among the remaining groups
  if (vr_remaining > 0) {
    
    # Determine which groups still have people waiting to be vaccinated
    groups_remaining <- which(vmax > 0)
    
    # Split remaining vaccinations proportionally among remaining groups
    priority_vaccinations[groups_remaining] <- priority_vaccinations[groups_remaining] + floor(
      pmin(vr_remaining * (vmax[groups_remaining]/sum(vmax[groups_remaining])),
           vmax[groups_remaining]
      )
    )
  }
  
  return(priority_vaccinations)
  
}
