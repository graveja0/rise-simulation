rm(list=ls())
source("./main_file.R")
source("./costs_simple.R")

## Look at summary statistics
results <- NULL
attributes <- NULL

#can overwrite default inputs
inputs$vHorizon <- 10
inputs$vN <- 100
inputs$vAge <- 40

#run and compile the simulation
for(preemptive in "None")
{
  for(reactive in c("None","Single"))
  {
    if(preemptive == "Single" && reactive == "Single") {next}
    inputs$vPreemptive <- preemptive
    inputs$vReactive   <- reactive
  cat("Running ", preemptive,"-", reactive,"\n")
  run <- exec.simulation(inputs) #obtain trajectory data
  run$preemptive <- preemptive
  run$reactive <- reactive
  at <- arrange(get_mon_attributes(env),name,key,time) #obtain attributes data
  at$preemptive <- preemptive
  at$reactive <- reactive

  if(is.null(results)) { results <- run } else  {results <- rbind(results, run)}
  if(is.null(attributes)) { attributes <- at } else  {attributes <- rbind(attributes, at)}
}}

#event counts
DT <- data.table(results)
summary <- DT[, .N, by = list(resource,preemptive,reactive)]
summary 

#cost and QALY
s1 <- cost.qaly(subset(results,preemptive=="None"&reactive=="None"),inputs) %>% mutate(strategy="None")
s2 <- cost.qaly(subset(results,preemptive=="None"&reactive=="Single"),inputs) %>% mutate(strategy="Reactive Single")

sum_costs <- rbind(s1,s2) %>% mutate(ICER = (dCOST-dCOST[1])/(dQALY-dQALY[1])) 
sum_costs
ic <- icer(sum_costs)
p <- ceplane(ic)