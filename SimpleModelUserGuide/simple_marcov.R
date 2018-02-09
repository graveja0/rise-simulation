# library(heemod)
# library(diagram)

rm(list=ls())
# fake risks for now
# to add discounting
# to add time dependency for secular death
# to add strategy

# states: 
# H
# Atrans: one year disu for A, lifetime drug cost
# BF: B free, lifetime drug cost
# BStrans: B survive, onetime cost for BS, lifetime disu for B, lifetime drug cost
# BS: B survive, lifetime disu for B, lifetime drug cost
# BDtrans: B death, onetime cost for BD
# D: both secular death and B death

#transition matrix
mat_stand <- define_transition(
        state_names = c("H","Atrans","BFree","BStrans","BS","BDtrans","D"),
        C,0.5,0,0,0,0,0.1,
        0,0,C,0.2,0,0.05,0.1,
        0,0,C,0.2,0,0.05,0.1,
        0,0,0,0,C,0,0.1,
        0,0,0,0,C,0,0.1,
        0,0,0,0,0,0,1,
        0,0,0,0,0,0,1
)

plot(mat_stand)

# param
costA = 10000
disuA = 0.05
costDrug = 365*0.5
costBS = 25000
disuB = 0.02
costBD = 15000


#states
state_H <- define_state(
        cost = 0,
        QALY = 1
)

state_Atrans <- define_state(
        cost = costA+costDrug,
        QALY = 1-disuA
)

state_BFree <- define_state(
        cost = costDrug,
        QALY = 1
)

state_BStrans <- define_state(
        cost = costBS+costDrug,
        QALY = 1-disuB
)

state_BS <- define_state(
        cost = costDrug,
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
strat_stand <- define_strategy(
        transition = mat_stand,
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
        strat_stand,
        cycles = 10,
        cost = cost,
        effect = QALY
)
