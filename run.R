rm(list=ls())

#can modify here
inputs.init <- list(
  vHorizon = 80,
  vN = 10000,
  vAge= 40,
  vN_PSA = 15
)

source("./sub-files/main_file.R")
source("./sub-files/costs_simple.R")

## Look at summary statistics
results <- NULL
attributes <- NULL

preemptive = "None"
reactive = "None"

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



# save(results,file='./run-data/result-toy.Rdata')
# save(at,file='./run-data/attributes-toy.Rdata')
# save(inputs,file='./run-data/inputs-toy.Rdata')
# Right now I am treating costs and disutilities as a constant, but in principle this could be done as a PSA too, by changing the first element
# from c(1) to seq(inputs.init$N_PSA).

#  s1 <- seq(inputs$vN_PSA) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"& reactive=="None" & aPSA_ID==.x),inputs=inputs)) %>% mutate(strategy="None")
#  s2 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None" & reactive=="Single" & aPSA_ID==.x),inputs=inputs)) %>% mutate(strategy="Reactive Single")
#  s3 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"&reactive=="Panel" & aPSA_ID==.x),inputs=inputs)) %>% mutate(strategy="Reactive Panel")
#  s4 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="Panel"&reactive=="None" & aPSA_ID==.x),inputs=inputs)) %>% mutate(strategy="Preemptive Panel")
# # 
#  sum_costs <- rbind(s1,s2,s3,s4) %>% arrange(dQALY) %>% mutate(ICER = (dCOST-dCOST[1])/(dQALY-dQALY[1])) 
#  sum_costs
