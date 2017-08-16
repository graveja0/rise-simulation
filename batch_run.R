rm(list=ls())
run.id <- "behavioral-global-test"
scenario = c("Panel","None")

preemptive = scenario[1]
reactive = scenario[2]

#can modify here
inputs.init <- list(
  vHorizon = 80,
  vN = 1500000,
  vAge= 40,
  vN_PSA = 200
)

source("./sub-files/main_file.R")
source("./sub-files/costs_simple.R")

## Look at summary statistics
results <- NULL
attributes <- NULL
select <- dplyr::select

inputs.init$vPreemptive <- preemptive
inputs.init$vReactive   <- reactive
cat("Running ", preemptive,"-", reactive,"\n")
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
save(results,file=file.path("./run-data/",paste0("results-",run.id,"-",preemptive,"-",reactive,".RData")))
save(attributes,file=file.path("./run-data/",paste0("attributes-",run.id,"-",preemptive,"-",reactive,".RData")))

DT <- results %>% arrange(aPSA_ID,name,start_time,end_time) %>% data.table()
summary <- DT[, .N, by = list(aPSA_ID,resource,preemptive,reactive)]
source("./sub-files/set-inputs.R")
save(drawn.parameter.values, file = file.path("./run-data/",paste0("parameter-values-",run.id,"-",preemptive,"-",reactive,".RData")))

DT <- results %>% arrange(aPSA_ID,name,start_time,end_time) %>% data.table()
summary <- DT[, .N, by = list(aPSA_ID,resource,preemptive,reactive)]
save(summary, file = file.path("./run-data/",paste0("summary-",run.id,"-",preemptive,"-",reactive,".RData")))

source("./sub-files/costs_simple.R")

# Get Overall Results
ce.results.overall <- cost.qaly(results,inputs) %>% mutate(strategy="None")
ce.results.psa <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results, aPSA_ID==.x),inputs=inputs,psa_id=.x)) 
ce.results <- list(overall=ce.results.overall,psa=ce.results.psa)
save(ce.results,file=file.path("./run-data/",paste0("ce-results-",run.id,"-",preemptive,"-",reactive,".RData")))


