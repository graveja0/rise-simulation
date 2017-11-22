days_till_A <- function(attrs, inputs)
{
  # Relative Risk
  rr <- 1
  
  # Baseline Risk
  days <- 365*inputs$vDurationA
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskA)*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["aTreat"]]==1) t2e <- 365*inputs$vHorizon + 1
  #if(t2e > days || attrs[["eventA"]] != 0) t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
} 


#for all genotyped patients through preemptive strategies (Panel or PREDICT), physician can choose to use or ignore the test results
#under reactive strategies, physician can also choose to order test for those not genotyped
A_reactive_strategy <- function(traj, inputs)
{
  if(inputs$vReactive == "None") 
  {
    traj #
  } else if (inputs$vReactive == "Single")
  {
    traj %>%
      branch(
        function(attrs) attrs[['aGenotyped']]+1,
        continue=c(TRUE, TRUE),
        trajectory("not have") %>%
          branch(
            function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
            continue=c(TRUE,TRUE),
            trajectory("not order") %>% timeout(0),
            trajectory("order reactive test") %>% set_attribute("aGenotyped", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test", 1)
            ), 
        trajectory("have test results") %>%  timeout(0)
      )
  } else if (inputs$vReactive == "Panel")
  {
    traj %>%
      branch(
        function(attrs) all_genotyped(attrs)+1,
        continue=c(TRUE, TRUE),
        trajectory("not panel tested") %>%
          branch(
            function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
            continue=c(TRUE,TRUE),
            trajectory("not order") %>% timeout(0),
            trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test", 1)
          ),           
        trajectory("panel tested") %>% timeout(0)
      )
  } else stop("Unhandled Reactive Strategy")
}

prescribe_drug <- function(traj,inputs) 
{
  traj %>%
  set_attribute("aDrug", function(attrs) 
    if(attrs[["aGene"]]==1 & attrs[["aGenotyped"]]==1 & 
       (attrs[['aOrdered_test']] == 1 | attrs[['aControlRead']]==1)) 
      {return(2)} else {return(1)}) %>%
    set_attribute("aTreat",1) %>%
    branch(
          function(attrs) attrs[["aDrug"]],
          continue=rep(TRUE,2),
          trajectory() %>% seize("rx"),
          trajectory() %>% seize("alt")
        )
}

event_A = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    set_attribute("aControlOrder",function() sample(0:1,1,prob=c(1- inputs$vProbabilityOrder,  inputs$vProbabilityOrder))) %>%
    set_attribute("aControlRead",function() sample(0:1,1,prob=c(1- inputs$vProbabilityRead,  inputs$vProbabilityRead))) %>% 
    A_reactive_strategy(inputs) %>%
    #assign drug
    prescribe_drug(inputs) %>%
    #event
    mark("A") %>% mark("A_c") %>%
    set_attribute("eventA",1) %>% #record occurance of A
    set_attribute("aRR_B",1) %>% #turn on B and adjust clock
    set_attribute("attB", function(attrs) now(env) + days_till_B(attrs,inputs)) 
}

days_till_B <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat"]]==1 & attrs[["aDrug"]]==2) inputs$vRR_B else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB)*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["eventA"]] != 1 || attrs[["eventB"]] != 0)
    t2e <- 365*inputs$vHorizon + 1

  return(t2e)
}

event_B = function(traj, inputs) 
{
  traj %>%
  mark("B") %>%
  set_attribute("eventB", 1) %>%
    branch(
      function(attrs) sample(1:2,1,prob=c(inputs$vFatalB,1-inputs$vFatalB)),
      continue = c(FALSE, TRUE),
      trajectory("Die")  %>% mark("B_Death") %>% terminate_simulation(),
      trajectory("Survive")  %>%  mark("B_Survive") 
    ) 
}



