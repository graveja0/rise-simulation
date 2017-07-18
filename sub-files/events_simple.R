days_till_A_a <- function(attrs, inputs)
{
  # Relative Risk
  rr <- 1
  
  # Baseline Risk
  days <- 365*inputs$vDurationA_a
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskA_a)*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["aTreat_a"]]==1) t2e <- 365*inputs$vHorizon + 1
  #if(t2e > days || attrs[["eventA"]] != 0) t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
} 

days_till_A_a_orig <- days_till_A_a

A_a_reactive_strategy <- function(traj, inputs)
{
  if(inputs$vReactive == "None") 
  {
    traj #
  } else if (inputs$vReactive == "Single")
  {
    traj %>%
      branch(
        function(attrs) attrs[['aGenotyped_a']]+1,
        continue=c(TRUE, TRUE),
        trajectory("not have") %>%
          branch(
            function(attrs) attrs[['aControlOrder_a']]+1, #use probability of ordering test
            continue=c(TRUE,TRUE),
            trajectory("not order") %>% timeout(0),
            trajectory("order reactive test") %>% set_attribute("aGenotyped_a", 1) %>% mark("single_test_a") %>% set_attribute("aOrdered_test_a", 1)
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
            function(attrs) attrs[['aControlOrder_a']]+1, #use probability of ordering test
            continue=c(TRUE,TRUE),
            trajectory("not order") %>% timeout(0),
            trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test_a", 1)
          ),           
        trajectory("panel tested") %>% timeout(0)
      )
  } else stop("Unhandled Reactive Strategy")
}

prescribe_drug_a <- function(traj,inputs) 
{
  traj %>%
    set_attribute("aDrug_a", function(attrs) 
      if(attrs[["aGene_a"]]==1 & attrs[["aGenotyped_a"]]==1 & 
         (attrs[['aOrdered_test_a']] == 1 | attrs[['aControlRead_a']]==1)) 
      {return(2)} else {return(1)}) %>%
    set_attribute("aTreat_a",1) %>%
    branch(
      function(attrs) attrs[["aDrug_a"]],
      continue=rep(TRUE,2),
      trajectory() %>% seize("rx_a"),
      trajectory() %>% seize("alt_a")
    )
}

event_A_a = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    set_attribute("aControlOrder_a",function() sample(0:1,1,prob=c(1- inputs$vProbabilityOrder_a,  inputs$vProbabilityOrder_a))) %>%
    set_attribute("aControlRead_a",function() sample(0:1,1,prob=c(1- inputs$vProbabilityRead_a,  inputs$vProbabilityRead_a))) %>% 
    A_a_reactive_strategy(inputs) %>%
    #assign drug
    prescribe_drug_a(inputs) %>%
    #event
    mark("A_a") %>% mark("A_c_a") %>%
    set_attribute("eventA_a",1) %>% #record occurance of A
    set_attribute("aRR_B_a",1) %>% #turn on B and adjust clock
    set_attribute("attB_a", function(attrs) now(env) + days_till_B_a(attrs,inputs)) 
}

days_till_B_a <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat_a"]]==1 & attrs[["aDrug_a"]]==2) inputs$vRR_B_a else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_a
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_a)*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(t2e > days || attrs[["eventA_a"]] != 1 || attrs[["eventB_a"]] != 0)
    t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
}

event_B_a = function(traj, inputs) 
{
  traj %>%
    mark("B_a") %>%
    set_attribute("eventB_a", 1) %>%
    branch(
      function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_a,1-inputs$vFatalB_a)),
      continue = c(FALSE, TRUE),
      trajectory("Die_a")  %>% mark("B_Death_a") %>% terminate_simulation(),
      trajectory("Survive_a")  %>%  mark("B_Survive_a") 
    ) 
}


sink("./temp/write-event-functions.R")
scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
                                  "days_till_A_TK <- function(attrs, inputs)
                                  {
                                  # Relative Risk
                                  rr <- 1
                                  
                                  # Baseline Risk
                                  days <- 365*inputs$vDurationA_TK
                                  
                                  # Convert To Probability 
                                  rate <- (- (log ( 1 - inputs$vRiskA_TK)*rr) / days)
                                  
                                  t2e <- rexp(1, rate)
                                  
                                  if(attrs[[\"aTreat_TK\"]]==1) t2e <- 365*inputs$vHorizon + 1
                                  
                                  return(t2e)
                                  }\n\n"
)) %>% unlist() %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
        "A_TK_reactive_strategy <- function(traj, inputs)
        {
          if(inputs$vReactive == \"None\") 
          {
          traj #
          } else if (inputs$vReactive == \"Single\")
          {
          traj %>%
          branch(
          function(attrs) attrs[[\'aGenotyped_TK\']]+1,
          continue=c(TRUE, TRUE),
          trajectory(\"not have\") %>%
          branch(
          function(attrs) attrs[[\'aControlOrder_TK\']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory(\"not order\") %>% timeout(0),
          trajectory(\"order reactive test\") %>% set_attribute(\"aGenotyped_TK\", 1) %>% mark(\"single_test_TK\") %>% set_attribute(\"aOrdered_test_TK\", 1)
          ), 
          trajectory(\"have test results\") %>%  timeout(0)
          )
          } else if (inputs$vReactive == \"Panel\")
          {
          traj %>%
          branch(
          function(attrs) all_genotyped(attrs)+1,
          continue=c(TRUE, TRUE),
          trajectory(\"not panel tested\") %>%
          branch(
          function(attrs) attrs[[\'aControlOrder_TK\']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory(\"not order\") %>% timeout(0),
          trajectory(\"order reactive test\") %>% panel_test() %>% set_attribute(\"aOrdered_test_TK\", 1)
          ),           
          trajectory(\"panel tested\") %>% timeout(0)
          )
          } else stop(\"Unhandled Reactive Strategy\")
          }\n\n"
)) %>% unlist() %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
"prescribe_drug_TK <- function(traj,inputs)
{
  traj %>%
    set_attribute(\"aDrug_TK\", function(attrs)
      if(attrs[[\"aGene_TK\"]]==1 & attrs[[\"aGenotyped_TK\"]]==1 &
                 (attrs[[\'aOrdered_test_TK\']] == 1 | attrs[[\'aControlRead_TK\']]==1))
                 {return(2)} else {return(1)}) %>%
                 set_attribute(\"aTreat_TK\",1) %>%
                 branch(
                 function(attrs) attrs[[\"aDrug_TK\"]],
                 continue=rep(TRUE,2),
                 trajectory() %>% seize(\"rx_TK\"),
                 trajectory() %>% seize(\"alt_TK\")
                 )
}\n\n"
)) %>% unlist() %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
"event_A_TK = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    set_attribute(\"aControlOrder_TK\",function() sample(0:1,1,prob=c(1- inputs$vProbabilityOrder_TK,  inputs$vProbabilityOrder_TK))) %>%
                   set_attribute(\"aControlRead_TK\",function() sample(0:1,1,prob=c(1- inputs$vProbabilityRead_TK,  inputs$vProbabilityRead_TK))) %>% 
                   A_TK_reactive_strategy(inputs) %>%
                   #assign drug
                   prescribe_drug_TK(inputs) %>%
                   #event
                   mark(\"A_TK\") %>% mark(\"A_c_TK\") %>%
                   set_attribute(\"eventA_TK\",1) %>% #record occurance of A
                   set_attribute(\"aRR_B_TK\",1) %>% #turn on B and adjust clock
                   set_attribute(\"attB_TK\", function(attrs) now(env) + days_till_B_TK(attrs,inputs)) 
}\n\n"
)) %>% unlist() %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
"days_till_B_TK <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[[\"aTreat_TK\"]]==1 & attrs[[\"aDrug_TK\"]]==2) inputs$vRR_B_TK else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_TK
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_TK)*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(t2e > days || attrs[[\"eventA_TK\"]] != 1 || attrs[[\"eventB_TK\"]] != 0)
  t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
  }\n\n
  "
  )) %>% unlist() %>% cat()

scenario.ids %>% purrr::map(~gsub("_TK",paste0("_",.x),
"event_B_TK = function(traj, inputs) 
{
  traj %>%
    mark(\"B_TK\") %>%
          set_attribute(\"eventB_TK\", 1) %>%
          branch(
          function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_TK,1-inputs$vFatalB_TK)),
          continue = c(FALSE, TRUE),
          trajectory(\"Die_TK\")  %>% mark(\"B_Death_TK\") %>% terminate_simulation(),
          trajectory(\"Survive_TK\")  %>%  mark(\"B_Survive_TK\") 
          ) 
}\n\n"
  )) %>% unlist() %>% cat()
          

sink()

source("./temp/write-event-functions.R")

#for all genotyped patients through preemptive strategies (Panel or PREDICT), physician can choose to use or ignore the test results
#under reactive strategies, physician can also choose to order test for those not genotyped



