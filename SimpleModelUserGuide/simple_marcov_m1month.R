library(heemod)
library(diagram)
library(tidyverse)

rm(list=ls())

# states: 
# H
# Atrans: one year disu for A, lifetime drug cost
# BF: B free, lifetime drug cost
# BStrans: B survive, onetime cost for BS, lifetime disu for B, lifetime drug cost
# BS: B survive, lifetime disu for B, lifetime drug cost
# BDtrans: B death, onetime cost for BD
# D: both secular death and B death

# parameters
load("raw_mortality.rda") #raw annual prob secular death, source data to fit death in simmer and numerical model.

doom <- function(value,cage) {
        ifelse(cage<max(lt$Age),value,0)
} # Once age reaches the maximum, pD becomes 1 and all other probs need to put 0.

doom2 <- function(value,cage) {
        ifelse(cage<=max(lt$Age),value,0)
} # Once age reaches the maximum, pD becomes 1 and all other probs need to put 0.

doom3 <- function(value) {
        ifelse(value>0,rescale_prob(p=value,to=1/12),0)
}

param <- define_parameters(
        costA = 10000,
        disuA = 0.05,
        costDrug = 365*0.5/12,
        costAlt = 365*5/12,
        costBS = 25000,
        disuB = 0.02,
        costBD = 15000,
        
        age_init = 40,
        age = age_init + markov_cycle,
        
        sd = map_dbl(age, function(x) doom2(lt$prob[lt$Age==x & lt$gender=="female"],x)),
        pD = map_dbl(sd, ~doom3(value=.x)),
        
        pA = map_dbl(age, function(x) doom(rescale_prob(p=0.1,from=10,to=1/12),x)),
        pB = map_dbl(age, function(x) doom(rescale_prob(p=0.02,from=1,to=1/12),x)),
        fatalB = map_dbl(age, function(x) doom(0.05,x)),
        pBS = pB*(1-fatalB),
        pBD = pB*fatalB,

        gene = 1, #0 or 1
        rr = 1-0.3*gene,
        
        #just for genotype strategy
        cDgenotype = map_dbl(gene, function(x) ifelse(x==0,costDrug,costAlt))
)

#transition matrix
mat_standard <- define_transition(
        state_names = c("H","A","BS","BD","D"),
        C,pA,0,0,pD,
        0,C,pBS,pBD,pD,
        0,0,C,0,pD,
        0,0,0,1,0,
        0,0,0,0,1
)

mat_genotype <- define_transition(
        state_names = c("H","A","BS","BD","D"),
        C,pA,0,0,pD,
        0,C,rr*pBS,rr*pBD,pD,
        0,0,C,0,pD,
        0,0,0,1,0,
        0,0,0,0,1
)

plot(mat_standard)

dr <- rescale_discount_rate(0.03,from=12,to=1) # discounting rate
#states
state_H <- define_state(
        cost = 0,
        QALY = discount(1,dr)
)

state_A <- define_state(
        cost = discount(dispatch_strategy(
                standard=costA*ifelse(state_time<=1,1,0)+costDrug,
                genotype=costA*ifelse(state_time<=1,1,0)+cDgenotype
        ), dr),
        QALY = discount(1-disuA*ifelse(state_time<=12,1,0),dr)
)


state_BS <- define_state(
        cost = discount(dispatch_strategy(
                standard=costBS*ifelse(state_time<=1,1,0)+costDrug,
                genotype=costBS*ifelse(state_time<=1,1,0)+cDgenotype
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_BD <- define_state(
        cost = discount(costBD*ifelse(state_time<=1,1,0),dr),
        QALY = 0
)

state_D <- define_state(
        cost = 0,
        QALY = 0
)

# binding 
strat_standard <- define_strategy(
        transition = mat_standard,
        H=state_H,
        A=state_A,
        BS=state_BS,
        BD=state_BD,
        D=state_D
)

strat_genotype <- define_strategy(
        transition = mat_genotype,
        H=state_H,
        A=state_A,
        BS=state_BS,
        BD=state_BD,
        D=state_D
)

# run
res_mod <- run_model(
        standard=strat_standard,
        genotype=strat_genotype,
        parameters = param,
        cycles = 24,
        cost = cost,
        effect = QALY,
        state_time_limit=12,
        method="beginning"
)


p <- data.frame(res_mod$eval_strategy_list$standard$parameters) #check parameters


### add gene prevalence
pop <- data.frame(
        gene=c(0,1),
        .weights=c(80,20)
)

res_h <- update(res_mod, newdata = pop)

# compare
# res_mod$run_model # old model
# res_h$updated_model # separate models and weights
res_h$combined_model$run_model # combined model
