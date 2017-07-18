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
           "skimr")
invisible(lapply(pkg, require, character.only = TRUE))


###initial inputs
epsilon <- 0.000000000001
inputs <- list(
  vAge = 40,
  vGender = 1,
  vGene_a = 0.2,
  
  vRiskA_a = 0.1,
  vDurationA_a = 10,
  
  vRiskB_a = 0.02,
  vDurationB_a = 1,
  vRR_B_a = 0.7,
  vFatalB_a = 0.05,
  
  vPreemptive = "None", # "None" or "Panel"
  vReactive = "None", # "None" or "Single" or "Panel"
  vHorizon  = 10,
  vN = 100,
  
  vProbabilityOrder_a = 0.5,
  vProbabilityRead_a = 0.5,
  
  disutilities = list(
    A_a = 0.05,
    B_Survive_a = 0.1,
    B_Death_a  = 1,
    secular_death = 1
  ),
  
  durations = list(
    A_a = 365
  ),
  
  type = list(
    A_a = 1,
    A_c_a = 0,
    B_Survive_a = 0,
    B_Death_a = 0,
    secular_death_a = 0
  ),
  
  costs = list(
    A_c_a = 10000, #separate event to capture A cost
    B_Survive_a = 25000, 
    B_Death_a = 15000,
    rx_a= 0.5,
    alt_a=5,
    single_test_a=100,
    panel_test=250
  )
)
inputs.orig <- inputs
source("./sub-files/read-in-scenarios.R")


###
###assign attributes
id <- 0
sink("./temp/initialize_patient_attributes.R")
paste("initialize_patient <- function(traj, inputs) 
      { 
      traj %>% 
      seize(\"time_in_model\") %>%
      set_attribute(\"aAgeInitial\", function() inputs$vAge) %>% 
      set_attribute(\"aAge\", function(attrs) attrs[['aAgeInitial']]) %>% 
      set_attribute(\"aGender\", function() inputs$vGender) %>% 
      ") %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aGene_",.x,"\", function() sample(1:2,1,prob=c(inputs$vGene_",.x,",1-inputs$vGene_",.x,")))"," %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aGenotyped_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"eventA_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"eventB_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aTreat_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aDrug_",.x,"\",1) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aOrdered_test_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aControlOrder_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aControlRead_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
paste("  set_attribute(\"last\" , 1)\n}") %>% cat()
sink()
source("./temp/initialize_patient_attributes.R")


#preemptive strategy
preemptive_strategy <- function(traj, inputs)
{
  #traj <- predict_draw(traj, inputs) # Always execute predict random draw to keep seeded random number
  # states the same
  
  # Note this doesn't have to use branch, because it's a global that every trajectory gets
  if        (inputs$vPreemptive == "None"     )
  {
    traj # Do nothing
  } else if (inputs$vPreemptive == "Panel"    )
  {
    traj %>% panel_test(inputs) 
  } else stop("Unhandled Preemptive Strategy")
}


########
# Define Panel Test attributes, functions

# JAG these will also need to be updated to see if you were genotyped for any of the scenarios.
all_genotyped <- function(attrs)
{
  attrs[['aGenotyped_a']]     == 1 #
    #attrs[['aGenotyped_CYP2C19']] == 1 &&  # Clopidogrel
    #attrs[['aGenotyped_Warfarin']] == 1    # Warfarin  
}

any_genotyped <- function(attrs)
{
  attrs[['aGenotyped_a']]     == 1 
  #attrs[['aGenotyped_CYP2C19']] == 1 ||
  #attrs[['aGenotyped_Warfarin']] == 1 
}

panel_test <- function(traj, inputs)
{
  traj %>% 
    set_attribute('aGenotyped_a', 1)  %>%
    #set_attribute('aGenotyped_CVD',     1)  %>%
    #set_attribute('aGenotyped_Warfarin', 1) %>%
    mark("panel_test")
}

sink("./temp/define-genotyping-attributes-functions.R")
paste0("all_genotyped <- function(attrs) {\n") %>% cat()
purrr::map(scenario.ids,~paste0("\t attrs[[\'aGenotyped_",.x,"\']] ==1")) %>% unlist() %>% paste(collapse=" && \n") %>% cat()
paste0("\n}\n") %>% cat()

paste0("any_genotyped <- function(attrs) {\n") %>% cat()
purrr::map(scenario.ids,~paste0("\t attrs[[\'aGenotyped_",.x,"\']] ==1")) %>% unlist() %>% paste(collapse=" || \n") %>% cat()
paste0("\n}\n") %>% cat()

paste0("panel_test <- function(traj,inputs) {\n") %>% cat()
paste0("\ttraj %>% \n") %>% cat()
purrr::map(scenario.ids,~paste0("\t set_attribute(\'aGenotyped_",.x,"\', 1)")) %>% unlist() %>% paste(collapse=" %>% \n") %>% cat()
paste0(" %>%\n\tmark(\"panel_test\")\n}\n") %>% cat()
sink()
source("./temp/define-genotyping-attributes-functions.R")

####
## Secular Death
source('./sub-files/event_secular_death.R')
source('./sub-files/events_simple.R')


####
## Cleanup 
cleanup_on_termination <- function(traj)
{
  traj %>% 
    release("time_in_model") %>%
    branch(
        function(attrs) attrs[["aTreat_a"]]+1,
        continue=rep(TRUE,2),
        trajectory() %>% timeout(0),
        trajectory() %>% 
          branch(
            function(attrs) attrs[["aDrug_a"]],
            continue=rep(TRUE,2),
            trajectory() %>% release("rx_a"),
            trajectory() %>% release("alt_a")
      )) 
}


sink("./temp/cleanup-function.R")
paste0("cleanup_on_termination <- function(traj)
       {
       traj %>% 
       ") %>% cat()
scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
                                  "branch(
                                  function(attrs) attrs[[\"aTreat_TK\"]]+1,
                                  continue=rep(TRUE,2),
                                  trajectory() %>% timeout(0),
                                  trajectory() %>% 
                                  branch(
                                  function(attrs) attrs[[\"aDrug_TK\"]],
                                  continue=rep(TRUE,2),
                                  trajectory() %>% release(\"rx_TK\"),
                                  trajectory() %>% release(\"alt_TK\")
                                  ))" )) %>% unlist() %>% paste0(" %>%") %>% cat()

paste0("\nrelease(\"time_in_model\")\n}") %>% cat()
sink()
source("./temp/cleanup-function.R")

terminate_simulation <- function(traj, inputs)
{
  traj %>%
    branch(
      function() 1, 
      continue=FALSE,
      trajectory() %>% cleanup_on_termination()
    )
}



###
#fill in event_registry

# the event registry will also need to be updated, much like we updated the inputs 

event_registry <- list(
  list(name          = "Secular Death",
       attr          = "aSecularDeathTime",
       time_to_event = days_till_death,
       func          = secular_death,
       reactive      = FALSE),
  list(name          = "Event A_a",
       attr          = "attA_a",
       time_to_event = days_till_A_a,
       func          = event_A_a,
       reactive      = FALSE),
  list(name          = "Event B_a",
       attr          = "attB_a",
       time_to_event = days_till_B_a,
       func          = event_B_a,
       reactive      = FALSE),
  

  
  
  
  list(name          = "Terminate at time horizonb",
       attr          = "aTerminate",
       time_to_event = function(attrs,inputs) 365.0*inputs$vHorizon,
       func          = terminate_simulation,
       reactive      = FALSE)
)


sink("./temp/create-event-registry.R")
paste0(
  "event_registry <- list(
  list(name          = \"Secular Death\", 
  attr          = \"aSecularDeathTime\", 
  time_to_event = days_till_death,  
  func          = secular_death,  
  reactive      = FALSE),  \n") %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
                                  "
                                  \t list(name          = \"Event A_TK\",
                                  \t \t attr          = \"attA_TK\",
                                  \t time_to_event = days_till_A_TK,
                                  \t \t func          = event_A_TK,
                                  \t \t reactive      = FALSE),
                                  \t list(name          = \"Event B_TK\",
                                  \t \t attr          = \"attB_TK\",
                                  \t \t time_to_event = days_till_B_TK,
                                  \t \t func          = event_B_TK,
                                  \t \t reactive      = FALSE),
                                  ")) %>% unlist() %>% cat()

paste0(
  "  list(name          = \"Terminate at time horizonb\",
  \t \t attr          = \"aTerminate\",
  \t \t time_to_event = function(attrs,inputs) 365.0*inputs$vHorizon,
  \t \t func          = terminate_simulation,
  \t \treactive      = FALSE)
) \n") %>% cat()

sink()
source("./temp/create-event-registry.R")


#### Counters


# JAG also needs to be updated
counters <- c(
  "time_in_model", 
  "secular_death",
  "A_a",
  "A_c_a", #separate event to capture A cost
  "B_a",
  "B_Survive_a",
  "B_Death_a",
  "treat_a",
  "rx_a",
  "alt_a",
  "single_test_a",
  
  "panel_test"
)
counters.scenarios <- scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),paste0(c("A_TK","A_c_TK","B_TK","B_Survive_TK","B_Death_TK","treat_TK","rx_TK","alt_TK","single_test_TK")))) %>% unlist()
counters <- c("time_in_model","secular_death",counters.scenarios,"panel_test")

source('./sub-files/event_main_loop_simple.R')


##########
# Start the clock!
exec.simulation <- function(inputs)
{
  set.seed(12345)
  env  <<- simmer("Simple")
  traj <- simulation(env, inputs)
  env %>% create_counters(counters)
  
  env %>%
    add_generator("patient", traj, at(rep(0, inputs$vN)), mon=2) %>%
    run(365*inputs$vHorizon+1) %>% # Simulate just past horizon
    wrap()
  
  get_mon_arrivals(env, per_resource = T)
}


