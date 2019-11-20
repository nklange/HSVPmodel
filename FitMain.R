library("tidyverse")
library("GA")
library("circular")
library("parallel")

load("Data/ol17_prepared.rda")

# Source functions ------------------

source("PrepDataFunctions.R") # prepares data for fitting: pulling out id set sizes etc into list
source("FittingAlgorithmsFunctions.R") # nlminb + GA wrapper
source("OptimisationFunctions.R") # objective functions
source("GASimulationFunctions.R") # VP by simulation (R + Rcpp)
source("nlminbNumIntFunctions.R")# VP by numerical integration (R + Rcpp)



# prepare data: nested by id for mclapply in mutate ----------------
# cvid, leftout to make it compatible with fitting cv data parallelized across training sets
datprep <- data_ol17 %>% 
  filter(exp == "ol17_e1") %>% 
  mutate(cvid = id,
         leftout = 0) %>% 
  group_nest(exp,cvid,keep=T)

# Set number of fitting runs ------------------------------

nrep <- 1

# Fit models --------------------------


test <- datprep[16,] %>%
  mutate(fit_ga_MK_RNplus = mclapply(data, fit_one_vp_ga, rep = nrep, objective = ll_vp_full,
                                model = "MK_RNplus", ll_fun = sim_fun_MK_RNplus,
                                lower = lowerRNplus, upper = upperRNplus, mc.cores = 1, mc.preschedule = FALSE) )  %>% 
  mutate(fit_nlminb_MK_RNplus = mclapply(data, fit_one_vp_nlminb, rep = nrep,
                                         startpar = fit_ga_MK_RNplus ,
                                         model = "MK_RNplus", objective = ll_vp_full, ll_fun = ll_vp3, type = "r",
                                         lower = .Machine$double.eps, upper = Inf, mc.cores = 1, mc.preschedule = FALSE))


# load some GA fitted data
load("limitsOl17.rda")

# extract parameters and test individual components
restest2$fit_ga_MK_RNplus[[16]] %>% 
  as.data.frame() %>% 
  .[1,1:4] %>% 
  unlist()


fit_one_vp_nlminb(data = datprep[16,]$data[[1]],model = "MK_RNplus", 
                  objective = ll_vp_full, ll_fun = ll_vp3, type = "r",
                  lower = .Machine$double.eps,rep = nrep,
                  startpar = restest2[16,]$fit_ga_MK_RNplus[[1]] %>% mutate(cvid = id,leftout=0))


#### test smaller parts
pars_tmp <- restest2[16,]$fit_ga_MK_RNplus[[1]] %>% 
  mutate(cvid = id,leftout=0) %>% 
  select(1:4) %>% 
  slice(1) %>% 
  unlist()

dp_tmp <- prep_data(datprep[16,]$data[[1]])

ll_vp_full(pars = pars_tmp,
           ll_fun = ll_vp3, type = "r",lower = 0,
           set_sizes = dp_tmp$set_sizes, 
           error_list = dp_tmp$datalist, model = "MK_RNplus")

# for setsize 1
ll_vp3(pars = pars, errors = errors, model = model,type=type)



int_fun_MK_RNplus(c(0,1e2,1e3,1e4,1e5,1e6),pars_tmp,0.5)

kappa <- 1e6

besselI(kappa,0,expon.scaled = TRUE) 
