###
###Define events (need to register event registry and counters)

days_till_A <- function(inputs)
{
        # Relative Risk
        rr <- 1
        
        # Baseline Risk
        days <- 365*inputs$vDurationA
        
        rate <- (- (log ( 1 - inputs$vRiskA)*rr) / days)
        
        t2e <- rexp(1, rate)
        
        if(get_attribute(env,"eventA") != 0) t2e <- 365*inputs$vHorizon + 1 #prevent reoccurence of event A
        
        return(t2e)
} 


#for all genotyped patients through preemptive strategies, physician can choose to use or ignore the test results
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
                                function() get_attribute(env,"aGenotyped")+1,
                                continue=c(TRUE, TRUE),
                                trajectory("not have") %>%
                                        branch(
                                                function() get_attribute(env,"aControlOrder")+1, #use probability of ordering test
                                                continue=c(TRUE,TRUE),
                                                trajectory("not order") %>% timeout(0), 
                                                trajectory("order reactive test") %>% 
                                                        set_attribute("aGenotyped", 1) %>% mark("single_test") %>% set_attribute("aOrdered_test", 1)
                                        ), 
                                trajectory("have test results") %>%  timeout(0)
                        )
        } else stop("Unhandled Reactive Strategy")
}

prescribe_drug <- function(traj,inputs) 
{
        traj %>%
                set_attribute("aDrug", function() 
                        if(get_attribute(env,"aGene")==1 & get_attribute(env,"aGenotyped")==1 & 
                           (get_attribute(env,"aOrdered_test") == 1 | get_attribute(env,"aControlRead")==1)) 
                        {return(2)} else {return(1)}) %>%
                branch(
                        function() get_attribute(env,"aDrug"),
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
                #different implementation of genetic testing based on strategies and physician behavioral parameters
                A_reactive_strategy(inputs) %>%
                #assign different treatment based on genetic type and whether tested or not 
                prescribe_drug(inputs) %>%
                #event
                mark("A") %>% mark("A_c") %>%
                set_attribute("eventA",1) %>% #record occurance of A
                #adjust clock for event B
                set_attribute("attB", function() now(env) + days_till_B(inputs)) 
}

days_till_B <- function(inputs) 
{
        # Relative Risk
        rr <- if(get_attribute(env,"eventA")==1 & get_attribute(env,"aDrug")==2) inputs$vRR_B else 1.0
        
        # Baseline Risk
        days <- 365*inputs$vDurationB
        
        rate <- (- (log ( 1 - inputs$vRiskB)*rr) / days)
        
        t2e <- rexp(1, rate)
        
        if(get_attribute(env,"eventA") != 1 || get_attribute(env,"eventB") != 0)
                t2e <- 365*inputs$vHorizon + 1
        
        return(t2e)
}

event_B = function(traj, inputs) 
{
        traj %>%
                mark("B") %>%
                set_attribute("eventB", 1) %>%
                branch(
                        function() sample(1:2,1,prob=c(inputs$vFatalB,1-inputs$vFatalB)),
                        continue = c(FALSE, TRUE),
                        trajectory("Die")  %>% mark("B_Death") %>% terminate_simulation(),
                        trajectory("Survive")  %>%  mark("B_Survive") 
                ) 
}



