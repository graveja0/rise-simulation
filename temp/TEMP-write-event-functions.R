days_till_A_SC_A <- function(attrs, inputs)
                                  {
                                  # Relative Risk
                                  rr <- 1
                                  
                                  # Baseline Risk
                                  days <- 365*inputs$vDurationA_SC_A[attrs[['aPSA_ID']]]
                                  
                                  # Convert To Probability 
                                  rate <- (- (log ( 1 - inputs$vRiskA_SC_A[attrs[['aPSA_ID']]])*rr) / days)
                                  
                                  t2e <- rexp(1, rate)
                                  
                                  if(attrs[["aTreat_SC_A"]]==1) t2e <- 365*inputs$vHorizon + 1
                                  
                                  return(t2e)
                                  }

 days_till_A_SC_B <- function(attrs, inputs)
                                  {
                                  # Relative Risk
                                  rr <- 1
                                  
                                  # Baseline Risk
                                  days <- 365*inputs$vDurationA_SC_B[attrs[['aPSA_ID']]]
                                  
                                  # Convert To Probability 
                                  rate <- (- (log ( 1 - inputs$vRiskA_SC_B[attrs[['aPSA_ID']]])*rr) / days)
                                  
                                  t2e <- rexp(1, rate)
                                  
                                  if(attrs[["aTreat_SC_B"]]==1) t2e <- 365*inputs$vHorizon + 1
                                  
                                  return(t2e)
                                  }

 days_till_A_SC_C <- function(attrs, inputs)
                                  {
                                  # Relative Risk
                                  rr <- 1
                                  
                                  # Baseline Risk
                                  days <- 365*inputs$vDurationA_SC_C[attrs[['aPSA_ID']]]
                                  
                                  # Convert To Probability 
                                  rate <- (- (log ( 1 - inputs$vRiskA_SC_C[attrs[['aPSA_ID']]])*rr) / days)
                                  
                                  t2e <- rexp(1, rate)
                                  
                                  if(attrs[["aTreat_SC_C"]]==1) t2e <- 365*inputs$vHorizon + 1
                                  
                                  return(t2e)
                                  }

 days_till_A_SC_D <- function(attrs, inputs)
                                  {
                                  # Relative Risk
                                  rr <- 1
                                  
                                  # Baseline Risk
                                  days <- 365*inputs$vDurationA_SC_D[attrs[['aPSA_ID']]]
                                  
                                  # Convert To Probability 
                                  rate <- (- (log ( 1 - inputs$vRiskA_SC_D[attrs[['aPSA_ID']]])*rr) / days)
                                  
                                  t2e <- rexp(1, rate)
                                  
                                  if(attrs[["aTreat_SC_D"]]==1) t2e <- 365*inputs$vHorizon + 1
                                  
                                  return(t2e)
                                  }

A_SC_A_reactive_strategy <- function(traj, inputs)
        {
          if(inputs$vReactive == "None") 
          {
          traj #
          } else if (inputs$vReactive == "Single")
          {
          traj %>%
          branch(
          function(attrs) attrs[['aGenotyped_SC_A']]+1,
          continue=c(TRUE, TRUE),
          trajectory("not have") %>%
          branch(
          function(attrs) attrs[['aControlOrder_SC_A']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% set_attribute("aGenotyped_SC_A", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test_SC_A", 1)
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
          function(attrs) attrs[['aControlOrder_SC_A']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test_SC_A", 1)
          ),           
          trajectory("panel tested") %>% timeout(0)
          )
          } else stop("Unhandled Reactive Strategy")
          }

 A_SC_B_reactive_strategy <- function(traj, inputs)
        {
          if(inputs$vReactive == "None") 
          {
          traj #
          } else if (inputs$vReactive == "Single")
          {
          traj %>%
          branch(
          function(attrs) attrs[['aGenotyped_SC_B']]+1,
          continue=c(TRUE, TRUE),
          trajectory("not have") %>%
          branch(
          function(attrs) attrs[['aControlOrder_SC_B']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% set_attribute("aGenotyped_SC_B", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test_SC_B", 1)
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
          function(attrs) attrs[['aControlOrder_SC_B']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test_SC_B", 1)
          ),           
          trajectory("panel tested") %>% timeout(0)
          )
          } else stop("Unhandled Reactive Strategy")
          }

 A_SC_C_reactive_strategy <- function(traj, inputs)
        {
          if(inputs$vReactive == "None") 
          {
          traj #
          } else if (inputs$vReactive == "Single")
          {
          traj %>%
          branch(
          function(attrs) attrs[['aGenotyped_SC_C']]+1,
          continue=c(TRUE, TRUE),
          trajectory("not have") %>%
          branch(
          function(attrs) attrs[['aControlOrder_SC_C']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% set_attribute("aGenotyped_SC_C", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test_SC_C", 1)
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
          function(attrs) attrs[['aControlOrder_SC_C']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test_SC_C", 1)
          ),           
          trajectory("panel tested") %>% timeout(0)
          )
          } else stop("Unhandled Reactive Strategy")
          }

 A_SC_D_reactive_strategy <- function(traj, inputs)
        {
          if(inputs$vReactive == "None") 
          {
          traj #
          } else if (inputs$vReactive == "Single")
          {
          traj %>%
          branch(
          function(attrs) attrs[['aGenotyped_SC_D']]+1,
          continue=c(TRUE, TRUE),
          trajectory("not have") %>%
          branch(
          function(attrs) attrs[['aControlOrder_SC_D']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% set_attribute("aGenotyped_SC_D", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test_SC_D", 1)
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
          function(attrs) attrs[['aControlOrder_SC_D']]+1, #use probability of ordering test
          #function(attrs) attrs[['aControlOrder']]+1, #use probability of ordering test
          continue=c(TRUE,TRUE),
          trajectory("not order") %>% timeout(0),
          trajectory("order reactive test") %>% panel_test() %>% set_attribute("aOrdered_test_SC_D", 1)
          ),           
          trajectory("panel tested") %>% timeout(0)
          )
          } else stop("Unhandled Reactive Strategy")
          }

prescribe_drug_SC_A <- function(traj,inputs)
{
  traj %>%
    set_attribute("aDrug_SC_A", function(attrs)
      if(attrs[["aGene_SC_A"]]==1 & attrs[["aGenotyped_SC_A"]]==1 &
                 (attrs[['aOrdered_test_SC_A']] == 1 | attrs[['aControlRead_SC_A']]==1))
                  #(attrs[['aOrdered_test_SC_A']] == 1 | attrs[['aControlRead']]==1))
                 {return(2)} else {return(1)}) %>%
                 set_attribute("aTreat_SC_A",1) %>%
                 branch(
                 function(attrs) attrs[["aDrug_SC_A"]],
                 continue=rep(TRUE,2),
                 trajectory() %>% seize("rx_SC_A"),
                 trajectory() %>% seize("alt_SC_A")
                 )
}

 prescribe_drug_SC_B <- function(traj,inputs)
{
  traj %>%
    set_attribute("aDrug_SC_B", function(attrs)
      if(attrs[["aGene_SC_B"]]==1 & attrs[["aGenotyped_SC_B"]]==1 &
                 (attrs[['aOrdered_test_SC_B']] == 1 | attrs[['aControlRead_SC_B']]==1))
                  #(attrs[['aOrdered_test_SC_B']] == 1 | attrs[['aControlRead']]==1))
                 {return(2)} else {return(1)}) %>%
                 set_attribute("aTreat_SC_B",1) %>%
                 branch(
                 function(attrs) attrs[["aDrug_SC_B"]],
                 continue=rep(TRUE,2),
                 trajectory() %>% seize("rx_SC_B"),
                 trajectory() %>% seize("alt_SC_B")
                 )
}

 prescribe_drug_SC_C <- function(traj,inputs)
{
  traj %>%
    set_attribute("aDrug_SC_C", function(attrs)
      if(attrs[["aGene_SC_C"]]==1 & attrs[["aGenotyped_SC_C"]]==1 &
                 (attrs[['aOrdered_test_SC_C']] == 1 | attrs[['aControlRead_SC_C']]==1))
                  #(attrs[['aOrdered_test_SC_C']] == 1 | attrs[['aControlRead']]==1))
                 {return(2)} else {return(1)}) %>%
                 set_attribute("aTreat_SC_C",1) %>%
                 branch(
                 function(attrs) attrs[["aDrug_SC_C"]],
                 continue=rep(TRUE,2),
                 trajectory() %>% seize("rx_SC_C"),
                 trajectory() %>% seize("alt_SC_C")
                 )
}

 prescribe_drug_SC_D <- function(traj,inputs)
{
  traj %>%
    set_attribute("aDrug_SC_D", function(attrs)
      if(attrs[["aGene_SC_D"]]==1 & attrs[["aGenotyped_SC_D"]]==1 &
                 (attrs[['aOrdered_test_SC_D']] == 1 | attrs[['aControlRead_SC_D']]==1))
                  #(attrs[['aOrdered_test_SC_D']] == 1 | attrs[['aControlRead']]==1))
                 {return(2)} else {return(1)}) %>%
                 set_attribute("aTreat_SC_D",1) %>%
                 branch(
                 function(attrs) attrs[["aDrug_SC_D"]],
                 continue=rep(TRUE,2),
                 trajectory() %>% seize("rx_SC_D"),
                 trajectory() %>% seize("alt_SC_D")
                 )
}

event_A_SC_A = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    #set_attribute("aProbabilityOrder_SC_A",function() c(inputs$vProbabilityOrder_SC_A[attrs[['aPSA_ID']]],1-inputs$vProbabilityOrder_SC_A[attrs[['aPSA_ID']]])) %>% 
    # This version allows for differential ordering & reading by scenario. 
    #set_attribute("aControlOrder_SC_A",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder_SC_A[attrs[['aPSA_ID']]],inputs$vProbabilityOrder_SC_A[attrs[['aPSA_ID']]]))) %>%
     #             set_attribute("aControlRead_SC_A",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead_SC_A[attrs[['aPSA_ID']]],inputs$vProbabilityRead_SC_A[attrs[['aPSA_ID']]]))) %>% 

                   A_SC_A_reactive_strategy(inputs) %>%
                   #assign drug
                   prescribe_drug_SC_A(inputs) %>%
                   #event
                   mark("A_SC_A") %>% mark("A_c_SC_A") %>%
                   set_attribute("eventA_SC_A",1) %>% #record occurance of A
                   set_attribute("aRR_B_SC_A",1) %>% #turn on B and adjust clock
                   set_attribute("attB_SC_A", function(attrs) now(env) + days_till_B_SC_A(attrs,inputs)) 
}

 event_A_SC_B = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    #set_attribute("aProbabilityOrder_SC_B",function() c(inputs$vProbabilityOrder_SC_B[attrs[['aPSA_ID']]],1-inputs$vProbabilityOrder_SC_B[attrs[['aPSA_ID']]])) %>% 
    # This version allows for differential ordering & reading by scenario. 
    #set_attribute("aControlOrder_SC_B",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder_SC_B[attrs[['aPSA_ID']]],inputs$vProbabilityOrder_SC_B[attrs[['aPSA_ID']]]))) %>%
     #             set_attribute("aControlRead_SC_B",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead_SC_B[attrs[['aPSA_ID']]],inputs$vProbabilityRead_SC_B[attrs[['aPSA_ID']]]))) %>% 

                   A_SC_B_reactive_strategy(inputs) %>%
                   #assign drug
                   prescribe_drug_SC_B(inputs) %>%
                   #event
                   mark("A_SC_B") %>% mark("A_c_SC_B") %>%
                   set_attribute("eventA_SC_B",1) %>% #record occurance of A
                   set_attribute("aRR_B_SC_B",1) %>% #turn on B and adjust clock
                   set_attribute("attB_SC_B", function(attrs) now(env) + days_till_B_SC_B(attrs,inputs)) 
}

 event_A_SC_C = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    #set_attribute("aProbabilityOrder_SC_C",function() c(inputs$vProbabilityOrder_SC_C[attrs[['aPSA_ID']]],1-inputs$vProbabilityOrder_SC_C[attrs[['aPSA_ID']]])) %>% 
    # This version allows for differential ordering & reading by scenario. 
    #set_attribute("aControlOrder_SC_C",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder_SC_C[attrs[['aPSA_ID']]],inputs$vProbabilityOrder_SC_C[attrs[['aPSA_ID']]]))) %>%
     #             set_attribute("aControlRead_SC_C",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead_SC_C[attrs[['aPSA_ID']]],inputs$vProbabilityRead_SC_C[attrs[['aPSA_ID']]]))) %>% 

                   A_SC_C_reactive_strategy(inputs) %>%
                   #assign drug
                   prescribe_drug_SC_C(inputs) %>%
                   #event
                   mark("A_SC_C") %>% mark("A_c_SC_C") %>%
                   set_attribute("eventA_SC_C",1) %>% #record occurance of A
                   set_attribute("aRR_B_SC_C",1) %>% #turn on B and adjust clock
                   set_attribute("attB_SC_C", function(attrs) now(env) + days_till_B_SC_C(attrs,inputs)) 
}

 event_A_SC_D = function(traj, inputs) 
{
  traj %>% 
    #physician behavior
    #set_attribute("aProbabilityOrder_SC_D",function() c(inputs$vProbabilityOrder_SC_D[attrs[['aPSA_ID']]],1-inputs$vProbabilityOrder_SC_D[attrs[['aPSA_ID']]])) %>% 
    # This version allows for differential ordering & reading by scenario. 
    #set_attribute("aControlOrder_SC_D",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder_SC_D[attrs[['aPSA_ID']]],inputs$vProbabilityOrder_SC_D[attrs[['aPSA_ID']]]))) %>%
     #             set_attribute("aControlRead_SC_D",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead_SC_D[attrs[['aPSA_ID']]],inputs$vProbabilityRead_SC_D[attrs[['aPSA_ID']]]))) %>% 

                   A_SC_D_reactive_strategy(inputs) %>%
                   #assign drug
                   prescribe_drug_SC_D(inputs) %>%
                   #event
                   mark("A_SC_D") %>% mark("A_c_SC_D") %>%
                   set_attribute("eventA_SC_D",1) %>% #record occurance of A
                   set_attribute("aRR_B_SC_D",1) %>% #turn on B and adjust clock
                   set_attribute("attB_SC_D", function(attrs) now(env) + days_till_B_SC_D(attrs,inputs)) 
}

days_till_B_SC_A <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat_SC_A"]]==1 & attrs[["aDrug_SC_A"]]==2) inputs$vRR_B_SC_A[attrs[['aPSA_ID']]] else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_SC_A[attrs[['aPSA_ID']]]
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_SC_A[attrs[['aPSA_ID']]])*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["eventA_SC_A"]] != 1 || attrs[["eventB_SC_A"]] != 0)
  t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
  }


   days_till_B_SC_B <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat_SC_B"]]==1 & attrs[["aDrug_SC_B"]]==2) inputs$vRR_B_SC_B[attrs[['aPSA_ID']]] else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_SC_B[attrs[['aPSA_ID']]]
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_SC_B[attrs[['aPSA_ID']]])*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["eventA_SC_B"]] != 1 || attrs[["eventB_SC_B"]] != 0)
  t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
  }


   days_till_B_SC_C <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat_SC_C"]]==1 & attrs[["aDrug_SC_C"]]==2) inputs$vRR_B_SC_C[attrs[['aPSA_ID']]] else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_SC_C[attrs[['aPSA_ID']]]
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_SC_C[attrs[['aPSA_ID']]])*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["eventA_SC_C"]] != 1 || attrs[["eventB_SC_C"]] != 0)
  t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
  }


   days_till_B_SC_D <- function(attrs,inputs) 
{
  # Relative Risk
  rr <- if(attrs[["aTreat_SC_D"]]==1 & attrs[["aDrug_SC_D"]]==2) inputs$vRR_B_SC_D[attrs[['aPSA_ID']]] else 1.0
  
  # Baseline Risk
  days <- 365*inputs$vDurationB_SC_D[attrs[['aPSA_ID']]]
  
  # Convert To Probability 
  rate <- (- (log ( 1 - inputs$vRiskB_SC_D[attrs[['aPSA_ID']]])*rr) / days)
  
  t2e <- rexp(1, rate)
  
  if(attrs[["eventA_SC_D"]] != 1 || attrs[["eventB_SC_D"]] != 0)
  t2e <- 365*inputs$vHorizon + 1
  
  return(t2e)
  }


  event_B_SC_A = function(traj, inputs) 
{
  traj %>%
    mark("B_SC_A") %>%
          set_attribute("eventB_SC_A", 1) %>%
          branch(
          function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_SC_A[attrs[['aPSA_ID']]],1-inputs$vFatalB_SC_A[attrs[['aPSA_ID']]])),
          continue = c(FALSE, TRUE),
          trajectory("Die_SC_A")  %>% mark("B_Death_SC_A") %>% terminate_simulation(),
          trajectory("Survive_SC_A")  %>%  mark("B_Survive_SC_A") 
          ) 
}

 event_B_SC_B = function(traj, inputs) 
{
  traj %>%
    mark("B_SC_B") %>%
          set_attribute("eventB_SC_B", 1) %>%
          branch(
          function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_SC_B[attrs[['aPSA_ID']]],1-inputs$vFatalB_SC_B[attrs[['aPSA_ID']]])),
          continue = c(FALSE, TRUE),
          trajectory("Die_SC_B")  %>% mark("B_Death_SC_B") %>% terminate_simulation(),
          trajectory("Survive_SC_B")  %>%  mark("B_Survive_SC_B") 
          ) 
}

 event_B_SC_C = function(traj, inputs) 
{
  traj %>%
    mark("B_SC_C") %>%
          set_attribute("eventB_SC_C", 1) %>%
          branch(
          function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_SC_C[attrs[['aPSA_ID']]],1-inputs$vFatalB_SC_C[attrs[['aPSA_ID']]])),
          continue = c(FALSE, TRUE),
          trajectory("Die_SC_C")  %>% mark("B_Death_SC_C") %>% terminate_simulation(),
          trajectory("Survive_SC_C")  %>%  mark("B_Survive_SC_C") 
          ) 
}

 event_B_SC_D = function(traj, inputs) 
{
  traj %>%
    mark("B_SC_D") %>%
          set_attribute("eventB_SC_D", 1) %>%
          branch(
          function(attrs) sample(1:2,1,prob=c(inputs$vFatalB_SC_D[attrs[['aPSA_ID']]],1-inputs$vFatalB_SC_D[attrs[['aPSA_ID']]])),
          continue = c(FALSE, TRUE),
          trajectory("Die_SC_D")  %>% mark("B_Death_SC_D") %>% terminate_simulation(),
          trajectory("Survive_SC_D")  %>%  mark("B_Survive_SC_D") 
          ) 
}

