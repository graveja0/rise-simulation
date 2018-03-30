library(heemod)
library(diagram)
library(purrr)

rm(list=ls())
# to add discounting

# states (t means temporary): 
# 1.H (Healthy)

# A just happens:
# 2.AtBF (no B): disu for A, cost for A + drug cost
# 3.AtBSt (B survive): disu for A + B, cost for A+BS, drug cost
# 4.AtBDt (B die): death, cost for A+BD

# post-A
# 5.BF (no B) : drug cost
# 6.BSt (BS just happens): disu for B, cost for BS, drug cost
# 7.BS (post-A&B, with B): disu for B, drug cost
# 8.BDt (BD just happens): death, cost for BD

# 9.D: both secular death and B death

# parameters
load("raw_mortality.rda") #raw annual prob secular death, source data to fit death in simmer and numerical model.

doom <- function(value,cage) {
        ifelse(cage<max(lt$Age),value,0)
} # Once age reaches the maximum, pD becomes 1 and all other probs need to put 0.

doom2 <- function(value,cage) {
        ifelse(cage<=max(lt$Age),value,0)
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
        age = age_init + markov_cycle,
        
        pD = map_dbl(age, function(x) doom2(lt$prob[lt$Age==x & lt$gender=="female"],x)),
        
        pA = map_dbl(age, function(x) doom(rescale_prob(p=0.1,from=10),x)),
        pB = map_dbl(age, function(x) doom(0.02,x)),
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
        state_names = c("H","AtBF","AtBSt","AtBDt","BF","BSt","BS","BDt","D"),
        C,pA*(1-pB),pA*pBS,pA*pBD,0,0,0,0,pD,
        0,0,0,0,C,pBS,0,pBD,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,0,0,1,
        0,0,0,0,C,pBS,0,pBD,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,0,0,1,
        0,0,0,0,0,0,0,0,1
)




mat_genotype <- define_transition(
        state_names = c("H","AtBF","AtBSt","AtBDt","BF","BSt","BS","BDt","D"),
        C,pA*(1-rr*pB),pA*pBS*rr,pA*pBD*rr,0,0,0,0,pD,
        0,0,0,0,C,rr*pBS,0,rr*pBD,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,0,0,1,
        0,0,0,0,C,rr*pBS,0,rr*pBD,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,C,0,pD,
        0,0,0,0,0,0,0,0,1,
        0,0,0,0,0,0,0,0,1
)


plot(mat_standard)

dr <- 0.03 # discounting rate
#states
state_H <- define_state(
        cost = 0,
        QALY = discount(1,dr)
)

state_AtBF <- define_state(
        cost = discount(dispatch_strategy(
                standard=costA+costDrug,
                genotype=costA+cDgenotype
        ), dr),
        QALY = discount(1-disuA,dr)
)

state_AtBSt <- define_state(
        cost = discount(dispatch_strategy(
                standard=costA+costBS+costDrug,
                genotype=costA+costBS+cDgenotype
        ), dr),
        QALY = discount(1-disuA-disuB,dr)
)

state_AtBDt <- define_state(
        cost = discount(costA+costBD, dr),
        QALY = 0
)

state_BF <- define_state(
        cost = discount(dispatch_strategy(
                standard=costDrug,
                genotype=cDgenotype
        ), dr),
        QALY = discount(1,dr)
)

state_BSt <- define_state(
        cost = discount(dispatch_strategy(
                standard=costBS+costDrug,
                genotype=costBS+cDgenotype
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_BS <- define_state(
        cost = discount(dispatch_strategy(
                standard=costDrug,
                genotype=cDgenotype
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_BDt <- define_state(
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
        AtBF=state_AtBF,
        AtBSt=state_AtBSt,
        AtBDt=state_AtBDt,
        BF=state_BF,
        BSt=state_BSt,
        BS=state_BS,
        BDt=state_BDt,
        D=state_D
)

strat_genotype <- define_strategy(
        transition = mat_genotype,
        H=state_H,
        AtBF=state_AtBF,
        AtBSt=state_AtBSt,
        AtBDt=state_AtBDt,
        BF=state_BF,
        BSt=state_BSt,
        BS=state_BS,
        BDt=state_BDt,
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
