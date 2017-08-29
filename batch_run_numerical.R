rm(list=ls())
main.results <- FALSE
pkg = list("tidyverse","deSolve","ellipse","readxl","lhs")
invisible(lapply(pkg, require, character.only = TRUE))
run.id.base <- "vogi-nber"

mainDir <- "./run-data"
subDir <- run.id.base

if (!file.exists(file.path(mainDir,subDir))){
  dir.create(file.path(mainDir, subDir))
}

ss_death <- read.csv("./numerical/ss-death-2011.csv")
source("./numerical/model-n-simple.R")
source("./numerical/numerical-functions.R")
source("./sub-files/metamodeling-functions.R")
times    <- seq(0, 80, by = 2 / 365)

if (main.results) {
  #####################################################
  # Single-Gene Testing: Indication-Specific Scenarios
  #####################################################
  Get_Single_Result <- function(ws = "DL_AH_SL", scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx") 
  {
    config  <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=1)$parameter.draws
    nber.none <- model.run(config, 1, "none",times) %>% t() %>% data.frame()
    nber.single <- model.run(config, 1, "reactive-single",times) %>% t() %>% data.frame()
    out <- rbind(nber.none,nber.single) %>% 
      mutate(strategy= c("none","reactive.single")) %>% 
      select(strategy,dCOST,dQALY) %>% Get_ICER()
    return(out)
  }
  
  single.scenarios <- c("DL_AH_SL","DL_AL_SH","DH_AL_SL","DH_AH_SL","DL_AH_SH","DH_AL_SH","DH_AH_SH")
  single.results <- single.scenarios %>% purrr::map_df(~Get_Single_Result(.x),.id="ID")
  single.results$scenario = single.scenarios[as.numeric(single.results$ID)]
  single.results.pairwise <- single.results %>% filter(ICER>0) %>% select(scenario,ICER) %>% arrange(ICER) %>% mutate(ICER=prettyNum(round(ICER,0),big.mark=",",scientific=FALSE))
  
  ##################
  # Overall Results
  ##################
  ws = "nber-base-fixed"
  scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
  parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=1)
  config  <- parameter.values$parameter.draws
  str.none <- model.run(config, 1, "none",times) %>% t() %>% data.frame()
  str.single <- model.run(config, 1, "reactive-single",times) %>% t() %>% data.frame()
  str.panel <- model.run(config, 1, "reactive-panel",times) %>% t() %>% data.frame()
  str.prepanel <- model.run(config, 1, "preemptive-panel",times) %>% t() %>% data.frame()
  lifetime.results <- rbind(str.none,str.single,str.panel,str.prepanel) %>% 
    mutate(strategy= c("none","reactive.single","reactive.panel","preemptive.panel")) %>% 
    select(strategy,dCOST,dQALY) %>% Get_ICER()
  
  # 10-Year Horizon
  times10    <- seq(0, 10, by = 2 / 365)
  str10.none <- model.run(config, 1, "none",times10) %>% t() %>% data.frame()
  str10.single <- model.run(config, 1, "reactive-single",times10) %>% t() %>% data.frame()
  str10.panel <- model.run(config, 1, "reactive-panel",times10) %>% t() %>% data.frame()
  str10.prepanel <- model.run(config, 1, "preemptive-panel",times10) %>% t() %>% data.frame()
  tenyr.results <- rbind(str10.none,str10.single,str10.panel,str10.prepanel) %>% 
    mutate(strategy= c("none","reactive.single","reactive.panel","preemptive.panel")) %>% 
    select(strategy,dCOST,dQALY) %>% Get_ICER()
  }

########################
# Overall Results: PSA
########################

vN_PSA = 4

library(doParallel)
library(broom)
ws = "nber-psa"
if (!exists("ss")) ss <- "none"
scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=vN_PSA)
config  <- parameter.values$parameter.draws 

# None
registerDoParallel(cores = parallel::detectCores())
getDoParWorkers()
psa.run <- foreach(ii = seq(vN_PSA), .combine = rbind) %dopar%
{
  library(deSolve)
  library(dplyr)
  model.run(config, ii, ss,times)
}

psa <- add.params(run=psa.run,parameter.draws=parameter.values$parameter.draws,run.name=ss) %>%  rename(psa_id=iteration) %>% select(-possible,-fatal_b,-living)
save(psa,file=file.path(mainDir, subDir,paste0("nber-psa-",ss,".RData")))


