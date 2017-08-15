#1. Request an interactive job with salloc (see
#                                           http://www.accre.vanderbilt.edu/?page <- id=2154#salloc )
#                                           2. Start a screen session on the compute node
#                                           3. Run your R process within the screen
#                                           4. Close the screen session and log out of the compute node
#                                           5. To interact with your process again, ssh to the compute node and reattach
#                                           the screen.
# ssh vmp1018

rm(list=ls())

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

preemptive = "None"
reactive = "None"
select <- dplyr::select

for(preemptive in c("None","Panel"))
{
  #for(reactive in c("None","Single","Panel"))
  for(reactive in c("None","Panel"))
  {
    if(preemptive == "PREDICT" && reactive == "Panel") {next}
    if(preemptive == "PREDICT" && reactive == "Single") {next}
    if(preemptive == "Panel" && reactive == "Single") {next}
    if(preemptive == "Panel" && reactive == "Panel") {next}
    #if(preemptive == "None" && reactive == "None") {next}
    inputs.init$vPreemptive <- preemptive
    inputs.init$vReactive   <- reactive
  cat("Running ", preemptive,"-", reactive,"\n")
  #run <- seq(inputs.init$vN_PSA) %>% purrr::map_df(~exec.simulation(inputs.init,ii=.x),.id="PSA_ID")
  run <- exec.simulation(inputs.init)
  run$preemptive <- preemptive
  run$reactive <- reactive
  at <- arrange(get_mon_attributes(env),name,key,time)
  at$preemptive <- preemptive
  at$reactive <- reactive
  
  psa_id <- at %>% filter(key=="aPSA_ID") %>% select(name,aPSA_ID=value)
  run <- run %>% left_join(psa_id,"name")
  
  if(is.null(results)) { results <- run } else  {results <- rbind(results, run)}
  if(is.null(attributes)) { attributes <- at } else  {attributes <- rbind(attributes, at)}
}}

DT <- results %>% arrange(aPSA_ID,name,start_time,end_time) %>% data.table()
summary <- DT[, .N, by = list(aPSA_ID,resource,preemptive,reactive)]
source("./sub-files/set-inputs.R")
save(drawn.parameter.values, file = "./run-data/drawn-parameter-values.Rdata")

# Get Overall Results
s1 <- cost.qaly(subset(results,preemptive=="None"&reactive=="None"),inputs) %>% mutate(strategy="None")
s2 <- cost.qaly(subset(results,preemptive=="None"&reactive=="Single"),inputs) %>% mutate(strategy="Reactive Single")
s3 <- cost.qaly(subset(results,preemptive=="None"&reactive=="Panel"),inputs) %>% mutate(strategy="Reactive Panel")
s4 <- cost.qaly(subset(results,preemptive=="Panel"&reactive=="None"),inputs) %>% mutate(strategy="Preemptive Panel")

overall_summary <- rbind(s1,s2,s3,s4) %>% mutate(ICER = (dCOST-dCOST[1])/(dQALY-dQALY[1])) 
overall_summary
save(overall_summary,file="./run-data/overall-summary.Rdata")


# Get PSA REsults
s1.i <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"& reactive=="None" & aPSA_ID==.x),inputs=inputs,psa_id=.x)) %>%
  rename(dQALY_None = dQALY , dCOST_None = dCOST) 
s2.i <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="Panel"& reactive=="None" & aPSA_ID==.x),inputs=inputs,psa_id=.x)) %>% 
  rename(dQALY_Panel_None = dQALY , dCOST_Panel_None = dCOST) 
s3.i <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"& reactive=="Panel" & aPSA_ID==.x),inputs=inputs,psa_id=.x)) %>% rename(dQALY_None_Panel = dQALY , dCOST_None_Panel = dCOST) 
s4.i <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"& reactive=="Single" & aPSA_ID==.x),inputs=inputs,psa_id=.x)) %>% rename(dQALY_None_Single = dQALY , dCOST_None_Single = dCOST) 

Sim <- cbind(s1.i,s2.i,s3.i,s4.i,drawn.parameter.values$global)
save(Sim,file="./run-data/simulation-results.Rdata")


