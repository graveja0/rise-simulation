rm(list=ls())

#can modify here
inputs.init <- list(
  vHorizon = 80,
  vN = 10,
  vAge= 40,
  vN_PSA = 100
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
  for(reactive in c("None","Single","Panel"))
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

  if(is.null(results)) { results <- run } else  {results <- rbind(results, run)}
  if(is.null(attributes)) { attributes <- at } else  {attributes <- rbind(attributes, at)}
}}

DT <- data.table(results)
summary <- DT[, .N, by = list(resource,preemptive,reactive)]
source("./sub-files/set-inputs.R")

# Right now I am treating costs and disutilities as a constant, but in principle this could be done as a PSA too, by changing the first element
# from c(1) to seq(inputs.init$N_PSA).

s1 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"& reactive=="None"),inputs=inputs,ii=.x)) %>% mutate(strategy="None")
s2 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None" & reactive=="Single"),inputs=inputs,ii=.x)) %>% mutate(strategy="Reactive Single")
s3 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="None"&reactive=="Panel"),inputs=inputs,ii=.x)) %>% mutate(strategy="Reactive Panel")
s4 <- c(1) %>% purrr::map_df(~cost.qaly(subset(results,preemptive=="Panel"&reactive=="None"),inputs=inputs,ii=.x)) %>% mutate(strategy="Preemptive Panel")

sum_costs <- rbind(s1,s2,s3,s4) %>% arrange(dQALY) %>% mutate(ICER = (dCOST-dCOST[1])/(dQALY-dQALY[1])) 
sum_costs
