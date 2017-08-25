rm(list=ls())

pkg = list("simmer",
           "ggplot2",
           "reshape2",
           "plyr", #need to load this before "dplyr"
           "tidyr",
           "dplyr",
           "msm",
           "data.table",
           "deSolve",
           "scales",
           "knitr")

invisible(lapply(pkg, require, character.only = TRUE))

save.results = TRUE
batch.parameters = TRUE
Batches <- 10

run.id.base <- "vogi"
Scenarios <- list (c("None","Single"),c("None","None"),c("None","Panel"),c("Panel","None"))

mainDir <- "./run-data"
subDir <- run.id.base

if (!file.exists(file.path(mainDir,subDir))){
  dir.create(file.path(mainDir, subDir))
}

# Assign Inputs
for (ss in Scenarios)
{
  set.seed(123) # Need the same parameters for each scenario 
  for (batch in seq(Batches)) 
  {
  run.id <- paste0(run.id.base,"-",batch)
  scenario = c(ss[1],ss[2])
  scenario.file ="./simple-pgx-scenario-parameters-vogi.csv"
  
  preemptive = scenario[1]
  reactive = scenario[2]
  
  #can modify here
  inputs.init <- list(
    vHorizon = 80,
    vN = 500000,
    vAge= 40,
    vN_PSA = 200
  )
 
  source("./sub-files/read-in-scenarios.R")
  source("./sub-files/set-inputs.R")
  if (batch==1) cat(inputs$vProbabilityOrder_SC_A[1:10])
  cat("\n")
  assign(paste0("inputs",paste0(ss,collapse="."),batch),inputs)
  assign(paste0("drawn.parameter.values",paste0(ss,collapse="."),batch),drawn.parameter.values)
  }
}

for (ss in Scenarios)
{
  set.seed(123) # Need the same parameters for each scenario 
  for (batch in seq(Batches)) 
  {
    run.id <- paste0(run.id.base, "-", batch)
    scenario = c(ss[1], ss[2])
    
    
    preemptive = scenario[1]
    reactive = scenario[2]
    
  inputs <- get(paste0("inputs",paste0(ss,collapse="."),batch))
  drawn.parameter.values <- get(paste0("drawn.parameter.values",paste0(ss,collapse="."),batch))
  
  source("./sub-files/main_file.R"); 
  source("./sub-files/costs_simple.R")
  cat(inputs$vProbabilityOrder_SC_A)
  cat("\n\n")
  ## Look at summary statistics
  results <- NULL
  attributes <- NULL
  select <- dplyr::select
  
  inputs.init$vPreemptive <- preemptive
  inputs.init$vReactive   <- reactive
  cat("Running ", preemptive,"-", reactive,"-",batch,"\n")
  if (!save.results) cat("NOTE: NOT SAVING RESULTS\n")
  run <- exec.simulation(inputs.init)
  run$preemptive <- preemptive
  run$reactive <- reactive
  at <- arrange(get_mon_attributes(env),name,key,time)
  at$preemptive <- preemptive
  at$reactive <- reactive
  
  psa_id <- at %>% filter(key=="aPSA_ID") %>% select(name,aPSA_ID=value)
  run <- run %>% left_join(psa_id,"name")
  results <- run
  attributes <- at
  #save(results,file=file.path("./run-data/",paste0("results-",run.id,"-",preemptive,"-",reactive,".RData")))
  #save(attributes,file=file.path("./run-data/",paste0("attributes-",run.id,"-",preemptive,"-",reactive,".RData")))
  

  if (save.results) save(drawn.parameter.values, file = file.path(mainDir, subDir,paste0("parameter-values-",run.id,"-",preemptive,"-",reactive,".RData")))
  
  # DT <- results %>% arrange(aPSA_ID,name,start_time,end_time) %>% data.table()
  # summary <- DT[, .N, by = list(aPSA_ID,resource,preemptive,reactive)]
  # save(summary, file = file.path(mainDir, subDir,paste0("summary-",run.id,"-",preemptive,"-",reactive,".RData")))
  # 
  source("./sub-files/costs_simple.R")
  
  # Get Overall Results
  ce.results.overall <- cost.qaly(results,inputs) %>% mutate(strategy="None")
  ce.results.psa <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results, aPSA_ID==.x),inputs=inputs,psa_id=.x)) 
  ce.results <- list(overall=ce.results.overall,psa=ce.results.psa)
  if (save.results) save(ce.results,file=file.path(mainDir, subDir,paste0("ce-results-",run.id,"-",preemptive,"-",reactive,".RData")))
  
  }
}

