initialize_patient <- function(traj, inputs) 
      { 
      traj %>% 
      seize("time_in_model") %>%
      set_attribute("aAgeInitial", function() inputs$vAge) %>% 
      set_attribute("aAge", function(attrs) attrs[['aAgeInitial']]) %>% 
      set_attribute("aGender", function() inputs$vGender) %>% 
      set_attribute("aPSA_ID",function() sample(seq(inputs$vN_PSA),1)) %>% 
      set_attribute("aControlOrder",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityOrder[attrs[['aPSA_ID']]],inputs$vProbabilityOrder[attrs[['aPSA_ID']]]))) %>%
      set_attribute("aControlRead",function(attrs) sample(0:1,1,prob=c(1-inputs$vProbabilityRead[attrs[['aPSA_ID']]],inputs$vProbabilityRead[attrs[['aPSA_ID']]]))) %>% 
        set_attribute("aGene_SC_A", function(attrs) sample(1:2,1,prob=c(inputs$vGene_SC_A[attrs[['aPSA_ID']]],1-inputs$vGene_SC_A[attrs[['aPSA_ID']]]))) %>% 
   set_attribute("aGene_SC_B", function(attrs) sample(1:2,1,prob=c(inputs$vGene_SC_B[attrs[['aPSA_ID']]],1-inputs$vGene_SC_B[attrs[['aPSA_ID']]]))) %>% 
   set_attribute("aGene_SC_C", function(attrs) sample(1:2,1,prob=c(inputs$vGene_SC_C[attrs[['aPSA_ID']]],1-inputs$vGene_SC_C[attrs[['aPSA_ID']]]))) %>% 
   set_attribute("aGene_SC_D", function(attrs) sample(1:2,1,prob=c(inputs$vGene_SC_D[attrs[['aPSA_ID']]],1-inputs$vGene_SC_D[attrs[['aPSA_ID']]]))) %>% 
  set_attribute("aGenotyped_SC_A",0) %>% 
   set_attribute("aGenotyped_SC_B",0) %>% 
   set_attribute("aGenotyped_SC_C",0) %>% 
   set_attribute("aGenotyped_SC_D",0) %>% 
  set_attribute("eventA_SC_A",0) %>% 
   set_attribute("eventA_SC_B",0) %>% 
   set_attribute("eventA_SC_C",0) %>% 
   set_attribute("eventA_SC_D",0) %>% 
  set_attribute("eventB_SC_A",0) %>% 
   set_attribute("eventB_SC_B",0) %>% 
   set_attribute("eventB_SC_C",0) %>% 
   set_attribute("eventB_SC_D",0) %>% 
  set_attribute("aTreat_SC_A",0) %>% 
   set_attribute("aTreat_SC_B",0) %>% 
   set_attribute("aTreat_SC_C",0) %>% 
   set_attribute("aTreat_SC_D",0) %>% 
  set_attribute("aDrug_SC_A",1) %>% 
   set_attribute("aDrug_SC_B",1) %>% 
   set_attribute("aDrug_SC_C",1) %>% 
   set_attribute("aDrug_SC_D",1) %>% 
  set_attribute("aOrdered_test_SC_A",0) %>% 
   set_attribute("aOrdered_test_SC_B",0) %>% 
   set_attribute("aOrdered_test_SC_C",0) %>% 
   set_attribute("aOrdered_test_SC_D",0) %>% 
  set_attribute("aControlOrder_SC_A",0) %>% 
   set_attribute("aControlOrder_SC_B",0) %>% 
   set_attribute("aControlOrder_SC_C",0) %>% 
   set_attribute("aControlOrder_SC_D",0) %>% 
  set_attribute("aControlRead_SC_A",0) %>% 
   set_attribute("aControlRead_SC_B",0) %>% 
   set_attribute("aControlRead_SC_C",0) %>% 
   set_attribute("aControlRead_SC_D",0) %>% 
  set_attribute("last" , 1)
}