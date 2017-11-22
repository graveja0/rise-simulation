rm(list=ls())

pkg = list("simmer",
           "ggplot2",
           "reshape2",
           "plyr", #need to load this before "dplyr"
           "tidyr",
           "dplyr",
           "msm",
           "data.table",
           "deSolve")
invisible(lapply(pkg, require, character.only = TRUE))


###initial inputs
epsilon <- 0.000000000001
inputs <- list(
  vAge = 40,
  vGender = 1,
  vGene = 0.2,
  
  vRiskA = 0.1,
  vDurationA = 10,

  vRiskB = 0.02,
  vDurationB = 1,
  vRR_B = 0.7,
  vFatalB = 0.05,
  
  vPreemptive = "None", # "None" or "Panel"
  vReactive = "None", # "None" or "Single" or "Panel"
  vHorizon  = 10,
  vN = 100,
  
  vProbabilityOrder = 0.5,
  vProbabilityRead = 0.5,
  
  disutilities = list(
    A = 0.05,
    B_Survive = 0.1,
    B_Death  = 1,
    secular_death = 1
  ),
  
  durations = list(
    A = 365
  ),
  
  type = list(
    A = 1,
    A_c = 0,
    B_Survive = 0,
    B_Death = 0,
    secular_death = 0
  ),
  
  costs = list(
    A_c = 10000, #separate event to capture A cost
    B_Survive = 25000, 
    B_Death = 15000,
    rx= 0.5,
    alt=5,
    single_test=100,
    panel_test=250
  )
)



###
###assign attributes
id <- 0

initialize_patient <- function(traj, inputs)
{
  traj %>%
    seize("time_in_model") %>%
    #set_attribute("aID", function() { tmp <- id; id <<- id + 1; tmp }) %>%
    set_attribute("aAgeInitial", function() inputs$vAge) %>%
    set_attribute("aAge", function(attrs) attrs[['aAgeInitial']]) %>%
    set_attribute("aGender", function() inputs$vGender) %>%
    set_attribute("aGene", function() sample(1:2,1,prob=c(inputs$vGene,1-inputs$vGene))) %>% #1 - targeted, 2 - not
    set_attribute("aGenotyped", 0) %>% # 0 - not, 1 - yes
    set_attribute("eventA",0) %>%  # Event A 0=not experienced, 1=experienced
    set_attribute("eventB",0) %>%  # Event B 0=not experienced, 1=experienced
    set_attribute("aTreat",0) %>% # On Treatment 0 - N, 1 - Y
    set_attribute("aDrug",1) %>% # 1 - usual drug, 2 - alt 
    set_attribute("aOrdered_test", 0)  %>%    # Did a physician order a test this time 0 - no, 1 - yes
    set_attribute("aControlOrder", 0) %>% #control ordering test 0 - no, 1 - yes
    set_attribute("aControlRead", 0) #control reading test 0 - no, 1 - yes
}

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
all_genotyped <- function(attrs)
{
  attrs[['aGenotyped']]     == 1 #
    #attrs[['aGenotyped_CYP2C19']] == 1 &&  # Clopidogrel
    #attrs[['aGenotyped_Warfarin']] == 1    # Warfarin  
}

any_genotyped <- function(attrs)
{
  attrs[['aGenotyped']]     == 1 
    #attrs[['aGenotyped_CYP2C19']] == 1 ||
    #attrs[['aGenotyped_Warfarin']] == 1 
}

panel_test <- function(traj, inputs)
{
  traj %>% 
    set_attribute('aGenotyped', 1)  %>%
    #set_attribute('aGenotyped_CVD',     1)  %>%
    #set_attribute('aGenotyped_Warfarin', 1) %>%
    mark("panel_test")
}

####
## Secular Death
source('./event_secular_death.R')
source('./events_simple.R')


####
## Cleanup 
cleanup_on_termination <- function(traj)
{
  traj %>% 
    release("time_in_model") %>%
    branch(
        function(attrs) attrs[["aTreat"]]+1,
        continue=rep(TRUE,2),
        trajectory() %>% timeout(0),
        trajectory() %>% 
          branch(
            function(attrs) attrs[["aDrug"]],
            continue=rep(TRUE,2),
            trajectory() %>% release("rx"),
            trajectory() %>% release("alt")
      )) 
}

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
event_registry <- list(
  list(name          = "Secular Death",
       attr          = "aSecularDeathTime",
       time_to_event = days_till_death,
       func          = secular_death,
       reactive      = FALSE),
  list(name          = "Event A",
       attr          = "attA",
       time_to_event = days_till_A,
       func          = event_A,
       reactive      = FALSE),
  list(name          = "Event B",
       attr          = "attB",
       time_to_event = days_till_B,
       func          = event_B,
       reactive      = FALSE),
  list(name          = "Terminate at time horizonb",
       attr          = "aTerminate",
       time_to_event = function(attrs,inputs) 365.0*inputs$vHorizon,
       func          = terminate_simulation,
       reactive      = FALSE)
)

#### Counters
counters <- c(
  "time_in_model", 
  "A",
  "A_c", #separate event to capture A cost
  "B",
  "B_Survive",
  "B_Death",
  "treat",
  "secular_death",
  "rx",
  "alt",
  "single_test",
  "panel_test"
)

source('./event_main_loop_simple.R')


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


