rm(list=ls())

pkg = list("simmer",
           "reshape2",
           "plyr", 
           "tidyr",
           "dplyr",
           "msm",
           "data.table",
           "deSolve")
invisible(lapply(pkg, require, character.only = TRUE))

###
###Defalt inputs
source('./inputs.R')

###
###Assign initial attributes
epsilon <- 0.000000000001 #use for zero-risk


initialize_patient <- function(traj, inputs)
{
        traj %>%
                seize("time_in_model") %>%
                set_attribute("aAgeInitial", function() inputs$vAge) %>%
                set_attribute("aAge", function(attrs) attrs[['aAgeInitial']]) %>%
                set_attribute("aGender", function() inputs$vGender) %>%
                set_attribute("aGene", function() sample(1:2,1,prob=c(inputs$vGene,1-inputs$vGene))) %>% #1 - targeted, 2 - not
                set_attribute("aGenotyped", 0) %>% # 0 - not, 1 - yes
                set_attribute("eventA",0) %>%  # Event A 0=not experienced, 1=experienced
                set_attribute("eventB",0) %>%  # Event B 0=not experienced, 1=experienced
                set_attribute("aDrug",1) %>% # 1 - usual drug, 2 - alt 
                set_attribute("aOrdered_test", 0)  %>%    # Did a physician order a test this time 0 - no, 1 - yes
                set_attribute("aControlOrder", 0) %>% #control ordering test 0 - no, 1 - yes
                set_attribute("aControlRead", 0) #control reading test 0 - no, 1 - yes
}

###
###some global functions

#preemptive strategy (this function is called in the main loop file and execute right after initialize patient and events)
preemptive_strategy <- function(traj, inputs)
{
        # Note this doesn't have to use branch, because it's a global that every trajectory gets
        if        (inputs$vPreemptive == "None"     )
        {
                traj # Do nothing
        } else if (inputs$vPreemptive == "Single"    )
        {
                traj %>% single_test(inputs) 
        } else stop("Unhandled Preemptive Strategy")
}

single_test <- function(traj, inputs)
{
        traj %>% 
                set_attribute('aGenotyped', 1)  %>%
                mark("single_test")
}

# Cleanup 
cleanup_on_termination <- function(traj)
{
        traj %>% 
                release("time_in_model") %>% 
                branch(
                        function(attrs) attrs[["eventA"]]+1,
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
###Event functions
source('./event_secular_death.R')
source('./events_simple.R')

###
###Event registry
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

###
###Counters (cost and QALY related inputs need to match counter name)
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
        "single_test"
)


###
###Execution functions

source('./event_main_loop_simple.R')

exec.simulation <- function(inputs)
{
        set.seed(12345) ###random seed makes simulation replicable 
        env  <<- simmer("Simple")
        traj <- simulation(env, inputs)
        env %>% create_counters(counters)
        
        env %>%
                add_generator("patient", traj, at(rep(0, inputs$vN)), mon=2) %>%
                run(365*inputs$vHorizon+1) %>% # Simulate just past horizon
                wrap()
        
        get_mon_arrivals(env, per_resource = T)
}


