# library(heemod)
# library(diagram)
# library(purrr)

rm(list=ls())
# to add discounting

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

probSD <- function(x) {
        dt <- lt %>% filter(Age==x, gender=="female")
        return(dt$prob)
}

doom <- function(value,cage) {
        ifelse(cage<max(lt$Age),value,0)
} # Once age reaches the maximum, pD becomes 1 and all other probs need to put 0.

param <- define_parameters(
        costA = 10000,
        disuA = 0.05,
        costDrug = 365*0.5,
        costAlt = 365*5,
        costBS = 25000,
        disuB = 0.02,
        costBD = 15000,
        
        age_init = 40,
        age = age_init + markov_cycle - 1,
        
        pD = map_dbl(age, function(x) doom(lt$prob[lt$Age==x & lt$gender=="female"],x)),
        
        pA = map_dbl(age, function(x) doom(rescale_prob(p=0.1,from=10),x)),
        pB = map_dbl(age, function(x) doom(0.02,x)),
        fatalB = map_dbl(age, function(x) doom(0.05,x)),
        pBS = pB*(1-fatalB),
        pBD = pB*fatalB,

        gene = 1, #0 or 1
        rr = 1-0.3*gene     
)

#transition matrix
mat_standard <- define_transition(
        state_names = c("H","Atrans","BFree","BStrans","BS","BDtrans","D"),
        C,pA,0,0,0,0,pD,
        0,0,C,pBS,0,pBD,pD,
        0,0,C,pBS,0,pBD,pD,
        0,0,0,0,C,0,pD,
        0,0,0,0,C,0,pD,
        0,0,0,0,0,0,1,
        0,0,0,0,0,0,1
)

mat_genotype <- define_transition(
        state_names = c("H","Atrans","BFree","BStrans","BS","BDtrans","D"),
        C,pA,0,0,0,0,pD,
        0,0,C,rr*pBS,0,rr*pBD,pD,
        0,0,C,rr*pBS,0,rr*pBD,pD,
        0,0,0,0,C,0,pD,
        0,0,0,0,C,0,pD,
        0,0,0,0,0,0,1,
        0,0,0,0,0,0,1
)


plot(mat_standard)

dr <- 0.03 # discounting rate
#states
state_H <- define_state(
        cost = 0,
        QALY = discount(1,dr)
)

state_Atrans <- define_state(
        cost = discount(dispatch_strategy(
                standard=costA+costDrug,
                genotype=costA+costAlt
        ), dr),
        QALY = discount(1-disuA,dr)
)

state_BFree <- define_state(
        cost = discount(dispatch_strategy(
                standard=costDrug,
                genotype=costAlt
        ), dr),
        QALY = discount(1,dr)
)

state_BStrans <- define_state(
        cost = discount(dispatch_strategy(
                standard=costBS+costDrug,
                genotype=costBS+costAlt
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_BS <- define_state(
        cost = discount(dispatch_strategy(
                standard=costDrug,
                genotype=costAlt
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_BDtrans <- define_state(
        cost = discount(costBD,dr),
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
        Atrans=state_Atrans,
        BFree=state_BFree,
        BStrans=state_BStrans,
        BS=state_BS,
        BDtrans=state_BDtrans,
        D=state_D
)

strat_genotype <- define_strategy(
        transition = mat_genotype,
        H=state_H,
        Atrans=state_Atrans,
        BFree=state_BFree,
        BStrans=state_BStrans,
        BS=state_BS,
        BDtrans=state_BDtrans,
        D=state_D
)

# run
res_mod <- run_model(
        standard=strat_standard,
        genotype=strat_genotype,
        parameters = param,
        cycles = 100,
        cost = cost,
        effect = QALY
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
