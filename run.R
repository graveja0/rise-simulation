rm(list=ls())
source("./sub-files/main_file.R")
source("./sub-files/costs_simple.R")

## Look at summary statistics
results <- NULL
attributes <- NULL

#can modify here
inputs$vHorizon <- 10
inputs$vN <- 10
inputs$vAge <- 40

for(preemptive in c("None","Panel"))
{
  for(reactive in c("None","Single","Panel"))
  {
    if(preemptive == "PREDICT" && reactive == "Panel") {next}
    if(preemptive == "PREDICT" && reactive == "Single") {next}
    if(preemptive == "Panel" && reactive == "Single") {next}
    if(preemptive == "Panel" && reactive == "Panel") {next}
    #if(preemptive == "None" && reactive == "None") {next}
    inputs$vPreemptive <- preemptive
    inputs$vReactive   <- reactive
  cat("Running ", preemptive,"-", reactive,"\n")
  run <- exec.simulation(inputs)
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

s1 <- cost.qaly(subset(results,preemptive=="None"&reactive=="None"),inputs) %>% mutate(strategy="None")
s2 <- cost.qaly(subset(results,preemptive=="None"&reactive=="Single"),inputs) %>% mutate(strategy="Reactive Single")
s3 <- cost.qaly(subset(results,preemptive=="None"&reactive=="Panel"),inputs) %>% mutate(strategy="Reactive Panel")
s4 <- cost.qaly(subset(results,preemptive=="Panel"&reactive=="None"),inputs) %>% mutate(strategy="Preemptive Panel")

sum_costs <- rbind(s1,s2,s3,s4) %>% arrange(dQALY) %>% mutate(ICER = (dCOST-dCOST[1])/(dQALY-dQALY[1])) 
sum_costs
