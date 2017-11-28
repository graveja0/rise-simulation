###
###default inputs 

inputs <- list(
        
        ###1. parameters used in simulation
        #global control
        vHorizon  = 10, #how many years to simulate
        vN = 100, #how 
        
        #demographics (can be constant value or draw from distribution/dataset)
        vAge = 40,
        vGender = 1,
        vGene = 0.2,
        
        #event risks
        #Example: Event A happens at a 10% rate over 10 years
        vRiskA = 0.1,
        vDurationA = 10, 
        
        #Example: Event B happens at a 2% rate over 1 year with a case fatality of 5%
        vRiskB = 0.02,
        vDurationB = 1,
        vFatalB = 0.05,
        #Example: alternative drug users have a relative risk of 0.7
        vRR_B = 0.7,
        
        
        #strategies
        #Example: default setting is none-none (no genetic testing);
        #Other combinations include: none-single (single gene testing at the time of indication),
        #                            none-panel (panel gene testing at the time of indication),
        #                            panel-none (panel gene testing at the beginning of the simulation)
        vPreemptive = "None", # "None" or "Panel"
        vReactive = "None", # "None" or "Single" or "Panel"
        
        #physician bvehavior
        vProbabilityOrder = 0.5, #the probability for doctors to order and make use of a test at the time of indication
        vProbabilityRead = 0.5, #the probability for doctors to use information from existing test at the time of indication
        
        
        ###2.parameters used in post-simulation computation (event name needs to match counters)
        disutilities = list(
                A = 0.05, 
                B_Survive = 0.1,
                B_Death  = 1,
                secular_death = 1
        ),
        
        #here only list event that temporarily affects utility 
        durations = list(
                A = 365
        ),
        
        #differentiate events that permanently or temporarily affect utility
        #1 - temporary (need to specify in "durations" section above)
        #0 - permanent (do not record in "durations" section above)
        type = list(
                #Example: Event A incurs an one-time cost and temporary disutility,
                #so we set "A" as temporary to calcualte the disutility, 
                #and "A_c" as permanent to calculate the cost.
                #Need to record both events in "counters" session and event file. 
                A = 1,
                A_c = 0, 
                
                B_Survive = 0,
                B_Death = 0,
                secular_death = 0
        ),
        
        costs = list(
                A_c = 10000, #separate event to capture A cost, note that "A" is not present here to avoid double counting
                B_Survive = 25000, 
                B_Death = 15000,
                rx= 0.5, #daily cost  
                alt=5, #daily cost
                single_test=100
        )
        
)
