library(heemod)
library(diagram)
library(tidyverse)
library(readxl)

###############################
# try to replicate BRCA model

# convenience functions to compute transtion probs
# https://pierucci.org/heemod/reference/probability.html
###############################

### Typos
#1."Dead" column, ovarian cancer death used bc chart
#2."Cost" column, ooph cost did not multiply the cost

### System
#1. Confusing cost and QALY computation
#2. No transition from no intervention to intervention in Markov model

rm(list=ls())
################################
# test vs no test: change probs of getting intervention
# stats: none, mast, ooph, duet, breast, ovarian, postcancer, dead

# Gene type
# High Penetrance Finding
## 1 - Pathogenic/Likely Pathogenic
## 2 - Variant of Unknown Significance
## 3 - Non-Pathogenic
# Moderate Penetrance Finding
## 4 - Pathogenic/Likely Pathogenic
## 5 - Variant of Unknown Significance
## 6 - Non-Pathogenic
# Low Penetrance Finding
## 7 - Pathogenic/Likely Pathogenic
## 8 - Variant of Unknown Significance
## 9 - Non-Pathogenic
# 10 - No variant
gene_prev <- c(unlist(map(list(.0392,.0499,.0109),function(x) x*c(.09,.27,.63))),.901) #decision tree parameters table

bc_table <- read_excel(path = "RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CH5:DA20") %>% 
        select(1,2,5,14) %>% set_names("age","age_range","cancer","death")
oc_table <- read_excel(path = "RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CH23:DA38") %>% 
        select(1,2,5,14) %>% set_names("age","age_range","cancer","death")

life_table <- read_excel(path = "RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CD4:CE105") %>% set_names("age","death")

#age-based cancer rate and cancer mortality
p_cancer <- function(cancer,age,which) { #transition probabilities table 
        n <- ceiling((age-19)/5)+1 #which line to read
        if(age>=85) n = 15
        
        if(cancer=="bc") {
                r <- c(bc_table$cancer[n],bc_table$death[n])
        } 
        if(cancer=="oc") {
                r <- c(oc_table$cancer[n],oc_table$death[n]*1.5)
        }
        if(age>=100) {
                r <- c(0,0)
        }
        return(r[which])
}

#gene-based cancer rr w/o intervention
# rr_cancer <- function(gene) { #parameters table
#         if(gene %in% c(3,6,9,10)) rr=1
#         if(gene %in% c(2,5,8)) rr=1.1
#         if(gene == 1) rr=15
#         if(gene == 4) rr=5
#         if(gene == 7) rr=1.5
#         return(rr)
# }

#secular death rate
secular <- function(age) {
        if(age>100) 1
        else life_table[life_table$age==age,][[2]]
}

param <- define_parameters(
        age_init = 45,
        age = age_init + model_time,
        
        #intervention rates depend on gene type
        int_plp = 0,#0.25,
        int_vus = 0,#0.05,

        #cancer base rates depend on age, no intervention ~ gene, with intervention ~ type of intervention
        pbc = map_dbl(age,~p_cancer(cancer="bc",age=.x,which=1)),
        poc = map_dbl(age,~p_cancer(cancer="oc",age=.x,which=1)),
        rr_bc_mast = 0.1,
        rr_bc_ooph = 0.51,
        rr_bc_duet = 0.05,
        rr_oc_mast = 1,
        rr_oc_ooph = 0.16,
        rr_oc_duet = 0.16,
        rr1 = 15,
        rr2 = 1.1,
        rr4 = 5,
        rr5 = 1.1,
        rr7 = 1.5,
        rr8 = 1.1,
        
        psd = map_dbl(age,~secular(.x)), #secular death
        pbcd = map_dbl(age,~p_cancer(cancer="bc",age=.x,which=2)), #breast cancer death
        pocd = map_dbl(age,~p_cancer(cancer="bc",age=.x,which=2)), #ovarian cancer death (typo in excel)
        # ppostd = map_dbl(psd,function(x) min(1,x)),
        
        unone = 0.95,
        uint = 0.95,
        ucancer = 0.8,
        upostc = 0.85,
        dmast = 0.05,
        dooph = 0.05,
        dduet = 0.1,
        
        cmast = 12596,
        cooph = 8144,
        cbc = 75873,
        coc = 127995,
        cpostc = 7738,
        cdeath = 65403
)

states <- c("mast","ooph","duet","hp_plp","hp_vus","mp_plp","mp_vus","lp_plp","lp_vus","np","breast","ovarian","postcancer","cancerd","dead")
mat_not <- define_transition(
        C,0,0,0,0,0,0,0,0,0,pbc*rr_bc_mast,poc*rr_oc_mast,0,0,psd,
        0,C,0,0,0,0,0,0,0,0,pbc*rr_bc_ooph,poc*rr_oc_ooph,0,0,psd,
        0,0,C,0,0,0,0,0,0,0,pbc*rr_bc_duet,poc*rr_oc_duet,0,0,psd,
        int_plp,int_plp,int_plp,C,0,0,0,0,0,0,pbc*rr1,poc*rr1,0,0,psd,
        int_vus,int_vus,int_vus,0,C,0,0,0,0,0,pbc*rr2,poc*rr2,0,0,psd,
        int_plp,int_plp,int_plp,0,0,C,0,0,0,0,pbc*rr4,poc*rr4,0,0,psd,
        int_vus,int_vus,int_vus,0,0,0,C,0,0,0,pbc*rr5,poc*rr5,0,0,psd,
        int_plp,int_plp,int_plp,0,0,0,0,C,0,0,pbc*rr7,poc*rr7,0,0,psd,
        int_vus,int_vus,int_vus,0,0,0,0,0,C,0,pbc*rr8,poc*rr8,0,0,psd,
        0,0,0,0,0,0,0,0,0,C,pbc,poc,0,0,psd,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pbcd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pocd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,0,psd,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        
        state_names = states
)


mat_not
plot(mat_not)

r <- 0.03
state_mast <- define_state(cost=0,QALY=discount(uint,r))
state_ooph <- define_state(cost=0,QALY=discount(uint,r))
state_duet <- define_state(cost=0,QALY=discount(uint,r))
# state_mast <- define_state(cost=as.integer(state_time==1)*cmast,QALY=uint-as.integer(state_time==1)*dmast)
# state_ooph <- define_state(cost=as.integer(state_time==1)*cooph,QALY=uint-as.integer(state_time==1)*dooph)
# state_duet <- define_state(cost=as.integer(state_time==1)*(cmast+cooph),QALY=uint-as.integer(state_time==1)*dduet)
state_hp_plp <- define_state(cost=discount(cmast*0.5,r),QALY=discount(unone-0.25*0.2,r))
state_hp_vus <- define_state(cost=discount(cmast*0.1,r),QALY=discount(unone-0.05*0.2,r))
state_mp_plp <- define_state(cost=discount(cmast*0.5,r),QALY=discount(unone-0.25*0.2,r))
state_mp_vus <- define_state(cost=discount(cmast*0.1,r),QALY=discount(unone-0.05*0.2,r))
state_lp_plp <- define_state(cost=discount(cmast*0.5,r),QALY=discount(unone-0.25*0.2,r))
state_lp_vus <- define_state(cost=discount(cmast*0.1,r),QALY=discount(unone-0.05*0.2,r))
state_np <- define_state(cost=0,QALY=discount(unone,r))
state_breast <- define_state(cost=discount(cbc,r),QALY=discount(ucancer,r))
state_ovarian <- define_state(cost=discount(coc,r),QALY=discount(ucancer,r))
state_postcancer <- define_state(cost=discount(cpostc,r),QALY=discount(upostc,r))
state_cancerd <- define_state(cost=discount(cdeath,r),QALY=0)
# state_dead <- define_state(cost=discount(as.integer(state_time==1)*cdeath,r),QALY=0)
state_dead <- define_state(cost=0,QALY=0)
strat <- define_strategy(
        transition = mat_not,
        mast=state_mast,
        ooph=state_ooph,
        duet=state_duet,
        hp_plp=state_lp_plp,
        hp_vus=state_lp_vus,
        mp_plp=state_mp_plp,
        mp_vus=state_mp_vus,
        lp_plp=state_lp_plp,
        lp_vus=state_lp_vus,
        np=state_np,
        
        breast=state_breast,
        ovarian=state_ovarian,
        postcancer=state_postcancer,
        cancerd=state_cancerd,
        dead=state_dead
)


#for init state
init_no <- define_init(
        mast=.01*1000L,
        ooph=.01*1000L,
        duet=.01*1000L,
        hp_plp=.97*1000L*gene_prev[1],
        hp_vus=.97*1000L*gene_prev[2],
        mp_plp=.97*1000L*gene_prev[4],
        mp_vus=.97*1000L*gene_prev[5],
        lp_plp=.97*1000L*gene_prev[7],
        lp_vus=.97*1000L*gene_prev[8],
        np=.97*1000L*sum(gene_prev[c(3,6,9,10)])
)

init_test <- define_init(
        mast=1000L*(.05*sum(gene_prev[c(1,4,7)])+.03*sum(gene_prev[c(2,5,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
        ooph=1000L*(.03*sum(gene_prev[c(1,4,7)])+.02*sum(gene_prev[c(2,5,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
        duet=1000L*(.02*sum(gene_prev[c(1,2,4,5,7,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
        hp_plp=.9*1000L*gene_prev[1],
        hp_vus=.93*1000L*gene_prev[2],
        mp_plp=.9*1000L*gene_prev[4],
        mp_vus=.93*1000L*gene_prev[5],
        lp_plp=.9*1000L*gene_prev[7],
        lp_vus=.93*1000L*gene_prev[8],
        np=.97*1000L*(gene_prev[3]+gene_prev[6]+gene_prev[9]+gene_prev[10])
)

res_mod_none <- run_model(
        strat = strat,
        parameters = param,
        cycles = 57,
        cost = cost,
        effect = QALY,
        state_time_limit = 1,
        init = init_no
)


res_mod_test <- run_model(
        strat = strat,
        parameters = param,
        cycles = 57,
        cost = cost,
        effect = QALY,
        state_time_limit = 1,
        init = init_test
)

# t <- res_mod$eval_strategy_list$standard$transition
# c1 <- data.frame(t[[1]])
# names(c1) <- attributes(t)$state_names
# row.names(c1) <- attributes(t)$state_names
# c2 <- data.frame(t[[2]])
# names(c2) <- attributes(t)$state_names
# row.names(c2) <- attributes(t)$state_names
# 
# see <- 

ct1 <- res_mod_none$eval_strategy_list$strat$counts
ct2 <- res_mod_test$eval_strategy_list$strat$counts



## add gene prevalence
# pop1 <- data.frame(
#         test=rep(0,10),
#         gene=1:10,
#         .weights=gene_prev
# )
# 
# res_h1 <- update(res_mod, newdata = pop1)
# res_h1$combined_model$run_model # combined model
# 
# pop2 <- data.frame(
#         test=rep(1,10),
#         gene=1:10,
#         .weights=gene_prev
# )
# 
# res_h2 <- update(res_mod, newdata = pop2)
# res_h2$combined_model$run_model # combined model

# compare
# res_mod$run_model # old model
# res_h$updated_model # separate models and weights
