# library(heemod)
# library(diagram)

rm(list=ls())
# to add discounting
# to add time dependency for secular death (pD)
# how to incorporate population heterogeneity (e.g. gene prevalence)
# https://pierucci.org/heemod/articles/g_heterogeneity.html

# states: 
# H
# Atrans: one year disu for A, lifetime drug cost
# BF: B free, lifetime drug cost
# BStrans: B survive, onetime cost for BS, lifetime disu for B, lifetime drug cost
# BS: B survive, lifetime disu for B, lifetime drug cost
# BDtrans: B death, onetime cost for BD
# D: both secular death and B death

# parameters
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
        rr = 0.7,
        pD = 0.01  #time-varying
        
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
