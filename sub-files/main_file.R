
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

# contains <- dplyr::contains()
# select <- dplyr::select()

###initial inputs
epsilon <- 0.000000000001

# Read in the various scenarios and store in the inputs list object.
source("./sub-files/read-in-scenarios.R")


###
###assign attributes
id <- 0
sink("./temp/TEMP-initialize-patient-attributes.R")
paste("initialize_patient <- function(traj, inputs) 
      { 
      traj %>% 
      seize(\"time_in_model\") %>%
      set_attribute(\"aAgeInitial\", function() inputs$vAge) %>% 
      set_attribute(\"aAge\", function(attrs) attrs[['aAgeInitial']]) %>% 
      set_attribute(\"aGender\", function() inputs$vGender) %>% 
      set_attribute(\"aPSA_ID\",function() sample(seq(inputs$vN_PSA),1)) %>% 
      #set_attribute(\"aControlOrder\",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder[attrs[[\'aPSA_ID\']]],inputs$vProbabilityOrder[attrs[[\'aPSA_ID\']]]))) %>%
      #set_attribute(\"aControlRead\",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead[attrs[[\'aPSA_ID\']]],inputs$vProbabilityRead[attrs[[\'aPSA_ID\']]]))) %>% 
      ") %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aGene_",.x,"\", function(attrs) sample(1:2,1,prob=c(inputs$vGene_",.x,"[attrs[[\'aPSA_ID\']]],1-inputs$vGene_",.x,"[attrs[[\'aPSA_ID\']]])))"," %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aGenotyped_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"eventA_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"eventB_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aTreat_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aDrug_",.x,"\",1) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~paste("  set_attribute(\"aOrdered_test_",.x,"\",0) %>% \n",sep="")) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
                                  paste("  set_attribute(\"aControlOrder_",.x,"\",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder_TK[attrs[[\'aPSA_ID\']]],inputs$vProbabilityOrder_TK[attrs[[\'aPSA_ID\']]]))) %>% \n",sep=""))
                            ) %>% unlist() %>% cat()
scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
                                  paste("  set_attribute(\"aControlRead_",.x,"\",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead_TK[attrs[[\'aPSA_ID\']]],inputs$vProbabilityRead_TK[attrs[[\'aPSA_ID\']]]))) %>% \n",sep="")) 
                            ) %>% unlist() %>% cat()
paste("  set_attribute(\"last\" , 1)\n}") %>% cat()
sink()
source("./temp/TEMP-initialize-patient-attributes.R")


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
sink("./temp/TEMP-define-genotyping-attributes-functions.R")
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
source("./temp/TEMP-define-genotyping-attributes-functions.R")

####
## Secular Death
source('./sub-files/event_secular_death.R')
source('./sub-files/events_simple.R')


####
## Cleanup 
sink("./temp/TEMP-cleanup-function.R")
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
source("./temp/TEMP-cleanup-function.R")

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

sink("./temp/TEMP-create-event-registry.R")
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
source("./temp/TEMP-create-event-registry.R")


#### Counters
counters.scenarios <- scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),paste0(c("A_TK","A_c_TK","B_TK","B_Survive_TK","B_Death_TK","treat_TK","rx_TK","alt_TK")))) %>% unlist()
counters <- c("time_in_model","secular_death",counters.scenarios,"panel_test","single_test")


source('./sub-files/event_main_loop_simple.R')


##########
# Start the clock!
exec.simulation <- function(inputs)
{
  set.seed(12345)
  
  # risks.as.list <- setNames(split(t(drawn.parameter.values[["risk"]][ii,]), seq(nrow(t(drawn.parameter.values[["risk"]][ii,])))), colnames(drawn.parameter.values[["risk"]]))
  # inputs.main <- append(list(
  #   vAge = 40,
  #   vGender = 1,
  #   vPreemptive = "None", # "None" or "Panel"
  #   vReactive = "None", # "None" or "Single" or "Panel"
  #   vHorizon  = 10,
  #   vN = 100
  # ),risks.as.list)
  # disutility.as.list <- setNames(split(t(drawn.parameter.values[["disutility"]][ii,]), seq(nrow(t(drawn.parameter.values[["disutility"]][ii,])))), colnames(drawn.parameter.values[["disutility"]]))
  # disutilities = append(list(
  #   secular_death = 1
  # ),disutility.as.list)
  # 
  # duration.as.list <- setNames(split(t(drawn.parameter.values[["duration"]][ii,]), seq(nrow(t(drawn.parameter.values[["duration"]][ii,])))), colnames(drawn.parameter.values[["duration"]]))
  # durations = append(list(
  # ),duration.as.list)
  # 
  # type.as.list <- setNames(split(t(drawn.parameter.values[["type"]][ii,]), seq(nrow(t(drawn.parameter.values[["type"]][ii,])))), colnames(drawn.parameter.values[["type"]]))
  # type = append(list(
  #   secular_death = 0
  # ),type.as.list)
  # 
  # 
  # cost.as.list <- setNames(split(t(drawn.parameter.values[["cost"]][ii,]), seq(nrow(t(drawn.parameter.values[["cost"]][ii,])))), colnames(drawn.parameter.values[["cost"]]))
  # costs = append(list(
  #   panel_test=250
  # ),cost.as.list)
  # 
  # inputs <- append(append(inputs.init,inputs.main),list(disutilities=disutilities,durations=durations,type=type,costs=costs))
  # 
  
  risks.as.list <- setNames(split(t(drawn.parameter.values[["risk"]]), seq(nrow(t(drawn.parameter.values[["risk"]])))), colnames(drawn.parameter.values[["risk"]]))
  global.as.list <- setNames(split(t(drawn.parameter.values[["global"]]), seq(nrow(t(drawn.parameter.values[["global"]])))), colnames(drawn.parameter.values[["global"]]))
  
  inputs.main <- append(list(
    vAge = 40,
    vGender = 1,
    vPreemptive = "None", # "None" or "Panel"
    vReactive = "None", # "None" or "Single" or "Panel"
    vHorizon  = 10,
    vN = 100
  ),risks.as.list)
  inputs.main <- append(inputs.main,global.as.list)
  disutility.as.list <- setNames(split(t(drawn.parameter.values[["disutility"]]), seq(nrow(t(drawn.parameter.values[["disutility"]])))), colnames(drawn.parameter.values[["disutility"]]))
  disutilities = append(list(
    secular_death = 1
  ),disutility.as.list)
  
  duration.as.list <- setNames(split(t(drawn.parameter.values[["duration"]]), seq(nrow(t(drawn.parameter.values[["duration"]])))), colnames(drawn.parameter.values[["duration"]]))
  durations = append(list(
  ),duration.as.list)
  
  type.as.list <- setNames(split(t(drawn.parameter.values[["type"]]), seq(nrow(t(drawn.parameter.values[["type"]])))), colnames(drawn.parameter.values[["type"]]))
  type = append(list(
    secular_death = 0
  ),type.as.list)
  
  
  cost.as.list <- setNames(split(t(drawn.parameter.values[["cost"]]), seq(nrow(t(drawn.parameter.values[["cost"]])))), colnames(drawn.parameter.values[["cost"]]))
  costs = append(list(
    panel_test=drawn.parameter.values$global$panel_test,
    single_test = drawn.parameter.values$global$single_test
  ),cost.as.list)
  
  inputs <- append(append(inputs.init,inputs.main),list(disutilities=disutilities,durations=durations,type=type,costs=costs))
  
  env  <<- simmer("Simple")
  traj <- simulation(env, inputs)
  env %>% create_counters(counters)
  
  env %>%
    add_generator("patient", traj, at(rep(0, inputs$vN)), mon=2) %>%
    simmer::run(365*inputs$vHorizon+1) %>% # Simulate just past horizon
    wrap()
  
  get_mon_arrivals(env, per_resource = T)
}




