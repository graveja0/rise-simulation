# library(heemod)
# library(diagram)

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

param <- define_parameters(
        costA = 10000,
        disuA = 0.05,
        costDrug = 365*0.5,
        costAlt = 365*5,
        costBS = 25000,
        disuB = 0.02,
        costBD = 15000,
        
        pA = rescale_prob(p=0.1,from=10),
        pB = 0.02,
        fatalB = 0.05,
        pBS = pB*(1-fatalB),
        pBD = pB*fatalB,
        
        age_init = 40,
        age = age_init + markov_cycle - 1,

        pD = map_dbl(age, function(x) lt$prob[lt$Age==x & lt$gender=="female"]),

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


#states
state_H <- define_state(
        cost = 0,
        QALY = 1
)

state_Atrans <- define_state(
        cost = dispatch_strategy(
                standard=costA+costDrug,
                genotype=costA+costAlt),
        QALY = 1-disuA
)

state_BFree <- define_state(
        cost = dispatch_strategy(
                standard=costDrug,
                genotype=costAlt),
        QALY = 1
)

state_BStrans <- define_state(
        cost = dispatch_strategy(
                standard=costBS+costDrug,
                genotype=costBS+costAlt),
        QALY = 1-disuB
)

state_BS <- define_state(
        cost = dispatch_strategy(
                standard=costDrug,
                genotype=costAlt),
        QALY = 1-disuB
)

state_BDtrans <- define_state(
        cost = costBD,
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
        cycles = 10,
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
res_mod$run_model
res_h$updated_model
res_h$combined_model$run_model
