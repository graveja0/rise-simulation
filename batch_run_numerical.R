rm(list=ls())
vN_PSA = 8
pkg = list("tidyverse","deSolve","ellipse","readxl","lhs")
invisible(lapply(pkg, require, character.only = TRUE))
run.id.base <- "vogi-numerical"
Scenarios <- list (c("None","Single"),c("None","None"),c("None","Panel"),c("Panel","None"))

mainDir <- "./run-data"
subDir <- run.id.base

if (!file.exists(file.path(mainDir,subDir))){
  dir.create(file.path(mainDir, subDir))
}

ss_death <- read.csv("./numerical/ss-death-2011.csv")
source("./numerical/model-n-simple.R")
source("./numerical/numerical-functions.R")
times    <- seq(0, 80, by = 2 / 365)

# # Value of Genomic Information (Single Gene Test)
# ws = "vogi-single-scenario"
# scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
# parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=1)
# config  <- parameter.values$parameter.draws 
# vogi.none1 <- model.run(config, 1, "none",times) %>% t() %>% data.frame()
# vogi.single1 <- model.run(config, 1, "reactive-single",times) %>% t() %>% data.frame()
# 
# ws = "vogi-panel20"
# scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
# parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=1)
# config  <- parameter.values$parameter.draws 
# vogi.none20 <- model.run(config, 1, "none",times) %>% t() %>% data.frame()
# vogi.single20 <- model.run(config, 1, "reactive-single",times) %>% t() %>% data.frame()
# 
# lambda <- seq(0,200000,1)
# NMB.none1 <- (vogi.none1$dQALY *lambda - vogi.none1$dCOST)
# NMB.single1 <- (vogi.single1$dQALY *lambda - vogi.single1$dCOST)
# 
# NMB.none20 <- (vogi.none20$dQALY *lambda - vogi.none20$dCOST)
# NMB.single20 <- (vogi.single20$dQALY *lambda - vogi.single20$dCOST)
# 
# VOGI20 <- NMB.single20-NMB.none20
# VOGI1 <- NMB.single1-NMB.none1
# lambda[VOGI20>0][1]
# lambda[VOGI1>0][1]
# 
# VOGI1[which(lambda==100000)]
# VOGI20[which(lambda==100000)]
# 
# 
# # Add in testing costs
# ws = "single-scenario"
# scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
# parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=1)
# config  <- parameter.values$parameter.draws 
# none1 <- model.run(config, 1, "none",times) %>% t() %>% data.frame()
# single1 <- model.run(config, 1, "reactive-single",times) %>% t() %>% data.frame()
# NMB.none1 <- (vogi.none1$dQALY *lambda - vogi.none1$dCOST)
# NMB.single1 <- (vogi.single1$dQALY *lambda - vogi.single1$dCOST)
# VOGI1 <- NMB.single1-NMB.none1
# lambda[VOGI1>0][1]
# VOGI1[which(lambda==100000)]


library(doParallel)
library(broom)
ws = "single-scenario-psa"
scenario.file = "./nber-aug2017/rise-nber-parameters-vogi-psa.xlsx"
parameter.values <- draw.parameter.values(scenario.file =scenario.file,ws=ws,PSA.N=vN_PSA)
config  <- parameter.values$parameter.draws 

# None
registerDoParallel(cores = parallel::detectCores())
getDoParWorkers()
vogi.none.psa.run <- foreach(ii = seq(vN_PSA), .combine = rbind) %dopar%
{
  library(deSolve)
  library(dplyr)
  model.run(config, ii, "none",times)
}
vogi.none.psa <- add.params(run=vogi.none.psa.run,parameter.draws=parameter.values$parameter.draws,run.name="none") %>%  rename(psa_id=iteration) %>% select(-possible,-fatal_b,-living)
save(vogi.none.psa,file=file.path(mainDir, subDir,"vogi-none-psa.RData"))


# Reactive Single
registerDoParallel(cores = parallel::detectCores())
getDoParWorkers()
vogi.single.psa.run <- foreach(ii = seq(vN_PSA), .combine = rbind) %dopar%
{
  library(deSolve)
  library(dplyr)
  model.run(config, ii, "reactive-single",times)
}
vogi.single.psa <- add.params(run=vogi.single.psa.run,parameter.draws=parameter.values$parameter.draws,run.name="none") %>%  rename(psa_id=iteration) %>% select(-possible,-fatal_b,-living)
save(vogi.single.psa,file=file.path(mainDir, subDir,"vogi-single-psa.RData"))

registerDoParallel(cores = parallel::detectCores())
getDoParWorkers()
vogi.panel.psa.run <- foreach(ii = seq(vN_PSA), .combine = rbind) %dopar%
{
  library(deSolve)
  library(dplyr)
  model.run(config, ii, "reactive-panel",times)
}
vogi.panel.psa <- add.params(run=vogi.panel.psa.run,parameter.draws=parameter.values$parameter.draws,run.name="none") %>%  rename(psa_id=iteration) %>% select(-possible,-fatal_b,-living)
save(vogi.panel.psa,file=file.path(mainDir, subDir,"vogi-panel-psa.RData"))


registerDoParallel(cores = parallel::detectCores())
getDoParWorkers()
vogi.prepanel.psa.run <- foreach(ii = seq(vN_PSA), .combine = rbind) %dopar%
{
  library(deSolve)
  library(dplyr)
  model.run(config, ii, "preemptive-panel",times)
}
vogi.prepanel.psa <- add.params(run=vogi.prepanel.psa.run,parameter.draws=parameter.values$parameter.draws,run.name="none") %>%  rename(psa_id=iteration) %>% select(-possible,-fatal_b,-living)
save(vogi.prepanel.psa,file=file.path(mainDir, subDir,"vogi-prepanel-psa.RData"))


# 
# 
# ```{r}
# 
# 
# varying1 <-
#   suppressWarnings(ce.results %>% map_dbl( ~ var(., na.rm = TRUE)))
# varying <- names(varying1[which(varying1 > 0)])
# id.vars <-
#   c("psa_id", "strategy", varying[-grep("psa_id|strategy|QALY|COST", varying)])
# wide.fmla <-
#   as.formula(paste0(paste0(c("psa_id", varying[-grep("psa_id|strategy|QALY|COST", varying)]), collapse =
#                              "+"), "~variable+strategy"))
# 
# ce.results <-
#   ce.results %>% reshape2::melt(id.vars = id.vars) %>% mutate(psa_id=as.character(psa_id)) %>% reshape2::dcast(wide.fmla)
# 
# 
# source("../sub-files/metamodeling-functions.R")
# 
# #Determine the number of strategies
# Strategies <- Strategies2 <- c("none","reactive.single","reactive.panel","preemptive.panel")
# ndep <- length(Strategies)
# 
# #Create vector of variable names
# Names <- ce.results %>% data.frame() %>% select(-psa_id) %>% colnames()
# Parms <- ce.results[, -grep("psa_id|Iteration|QALY|COST", colnames(ce.results))]   
# #Get parameter names
# paramNames <- colnames(Parms)
# indep <- ncol(Parms)
# Outcomes <- ce.results[, grep("psa_id|strategy|QALY|COST", colnames(ce.results))] 
# 
# lambda <- 100000
# 
# for (i in seq(length(Strategies))) ce.results[,paste0("NHB_",Strategies2[i])] <- ce.results %>% select(-dplyr::contains("NHB"),-dplyr::contains("NMB")) %>% select_if(grepl(Strategies2[i],names(.))) %>% mutate_at(vars(contains("COST")),funs(-1*./lambda)) %>% mutate(NHB=rowSums(.)) %>% pull(NHB)
# 
# for (i in seq(length(Strategies))) ce.results[,paste0("NMB_",Strategies2[i])] <- ce.results %>% select(-contains("NHB"),-dplyr::contains("NMB")) %>% select_if(grepl(Strategies2[i],names(.))) %>% mutate_at(vars(contains("COST")),funs(-1*.)) %>% mutate_at(vars(contains("QALY")),funs(lambda*.)) %>%  mutate(NMB=rowSums(.)) %>% pull(NMB)
# 
# ce.results <- ce.results %>% mutate(Iteration=row_number())
# 
# NHB <- ce.results %>% dplyr::select(dplyr::contains("NHB"))
# NMB <- ce.results %>% dplyr::select(dplyr::contains("NMB"))
# 
# ```
# 
# 
# 
# ```{r}
# nmb.gg <- reshape2::melt(NMB,  
#                          variable.name = "Strategy", 
#                          value.name = "NMB")
# colnames(nmb.gg) <- c("Strategy","NMB")
# 
# ## Plot NMB for different strategies
# require(ggplot2)
# require(scales)  # For dollar labels
# require(grid)
# number_ticks <- function(n) {function(limits) pretty(limits, n)} 
# # Faceted plot by Strategy
# ggplot(nmb.gg, aes(x = NMB/1000)) +
#   geom_histogram(aes(y =..density..), col="black", fill = "gray") +
#   geom_density(color = "red") +
#   facet_wrap(~ Strategy, scales = "free_y") +
#   xlab("Net Monetary Benefit (NMB) x10^3") +
#   scale_x_continuous(breaks = number_ticks(5), labels = dollar) + 
#   scale_y_continuous(breaks = number_ticks(5)) + 
#   theme_bw()
# 
# ```
# 
# 
# ```{r}
# n.sim <- dim(NMB)[1]
# 
# Strategy1 <- "none"
# Strategy2 <- "reactive.single"
# 
# inmb <- NMB %>% select_at(vars(contains(Strategy1),contains(Strategy2))) %>% mutate_at(vars(contains(Strategy2)),funs(-1*.)) %>%  mutate(Net_NMB=rowSums(.)) %>% 
#   mutate(Simulation = row_number()) %>% select(Simulation,Net_NMB)
# 
# #### Incremental NMB (INMB) ####
# # Calculate INMB of B vs A
# # Only B vs A but we could have plotted all combinations
# 
# ## Format data frame suitably for plotting
# inmb.gg <- reshape2::melt(inmb, id.vars = "Simulation", 
#                           variable.name = "Comparison", 
#                           value.name = "INMB")
# colnames(inmb.gg) <- c("Simulation","Comparison","INMB")
# txtsize<-16
# ## Plot INMB
# ggplot(inmb.gg, aes(x = INMB/1000)) +
#   geom_histogram(aes(y =..density..), col="black", fill = "gray") +
#   geom_density(color = "red") +
#   geom_vline(xintercept = 0, col = 4, size = 1.5, linetype = "dashed") +
#   facet_wrap(~ Comparison, scales = "free_y") +
#   xlab("Incremental Net Monetary Benefit (INMB) in thousand $") +
#   scale_x_continuous(breaks = number_ticks(5), limits = c(-100, 100)) + 
#   scale_y_continuous(breaks = number_ticks(5)) + 
#   theme_bw(base_size = 14)
# 
# ```
# 
# 
# ```{r}
# #### Loss Matrix ####  
# # Find optimal strategy (d*) based on the highest expected NMB
# d.star <- which.max(colMeans(NMB))
# d.star
# 
# # Or without iterating (much faster!)
# loss <- as.matrix(NMB - NMB[, d.star])
# head(loss)
# 
# #### EVPI ####
# ## Find maximum loss overall strategies at each state of the world 
# ## (i.e., PSA sample)
# require(matrixStats)
# max.loss.i <- rowMaxs(loss)
# head(max.loss.i)
# ## Average across all states of the world
# evpi <- mean(max.loss.i)
# evpi
# 
# 
# lambda_range<-c(1000,350000,50)
# ```
# 
# 
