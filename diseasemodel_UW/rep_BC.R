library(heemod)
library(diagram)
library(tidyverse)
library(readxl)

###############################
# try to replicate BRCA model

# convenience functions to compute transtion probs
# https://pierucci.org/heemod/reference/probability.html
###############################

### Compare
# event counts and QALY almost match (tiny issue with summing up to 1 making postcancer and death count smaller)
# costs cannot match (VLOOKUP oc death missing col number, cancer death cost not using the same term)

### Questions
#Under test strategy, actual transition from no intervention to intervention states; under none strategy, no actual transition but calculate in QALY
#Negative values (transition probs out of no-intervention >1 under test strategy) [cap in heemod]
#"Cost" column, missing ooph cost, cancer death cost inconsistent with counts, VLOOKUP oc death missing col number
# ("Dead" column  took last-year cancer counts and current-year mortality, "Cost" column took current-year cancer counts)
#"Dead" column, ovarian cancer death referenceing breast cancer chart


rm(list=ls())
################################
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

# 

################################

###0. Read in params from excel
params <- read_excel(path = "~/Box Sync/*temp/RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "B37:C145",col_types = c("text","numeric")) %>% set_names("param","value") %>%
        filter(!is.na(value))
bc_table <- read_excel(path = "/Box Sync/*temp/RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CH5:DA20") %>% 
        select(1,2,5,14) %>% set_names("age","age_range","cancer","death")
oc_table <- read_excel(path = "/Box Sync/*temp/RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CH23:DA38") %>% 
        select(1,2,5,14) %>% set_names("age","age_range","cancer","death")
life_table <- read_excel(path = "/Box Sync/*temp/RISE Model_BC.xlsm", sheet = "Breast_Ovarian", range = "CD4:CE105") %>% set_names("age","death")

###1. Define inputs
# gene type prevalence
gene_prev <- c(unlist(map(params$value[1:3]/10,function(x) x*params$value[4:6])),.9)

# age-based cancer rate and cancer mortality
p_cancer <- function(cancer,age,which) {
        n <- floor((age-19)/5)+1 #which line to read
        if(age>=85) n = 14
        
        if(cancer=="bc") {
                r <- c(bc_table$cancer[n],bc_table$death[n])
        } 
        if(cancer=="oc") {
                r <- c(oc_table$cancer[n],oc_table$death[n])
        }
        if(age>=100) {
                r <- c(0,0)
        }
        return(r[which])
}

#secular death rate
secular <- function(age) {
        if(age>100) 1
        else life_table[life_table$age==age,][[2]]
}

#cap psd
capf <- function(psd,pbc,poc,int,rr) {
        p <- min(psd,1-rr*(poc+pbc)-3*int)
        p
}

param <- define_parameters(
        age_init = 45,
        age = age_init + model_time,
        
        #intervention rates depend on gene type (only work under test strategy)
        int_plp = 0.25,
        int_vus = 0.05,
        
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
        pocd = map_dbl(age,~p_cancer(cancer="bc",age=.x,which=2)), #ovarian cancer death (typo in excel, should be "oc")
        ppostd = map_dbl(psd,function(x) min(1,x*1.5)),
        
        psd_cap_1 = pmap_dbl(list(psd,pbc,poc,int_plp),~capf(..1,..2,..3,..4,15)),
        psd_cap_2 = pmap_dbl(list(psd,pbc,poc,int_vus),~capf(..1,..2,..3,..4,1.1)),
        psd_cap_4 = pmap_dbl(list(psd,pbc,poc,int_plp),~capf(..1,..2,..3,..4,5)),
        psd_cap_5 = pmap_dbl(list(psd,pbc,poc,int_vus),~capf(..1,..2,..3,..4,1.1)),
        psd_cap_7 = pmap_dbl(list(psd,pbc,poc,int_plp),~capf(..1,..2,..3,..4,1.5)),
        psd_cap_8 = pmap_dbl(list(psd,pbc,poc,int_vus),~capf(..1,..2,..3,..4,1.1)),
        
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
mat_test <- define_transition(
        C,0,0,0,0,0,0,0,0,0,pbc*rr_bc_mast,poc*rr_oc_mast,0,0,psd,
        0,C,0,0,0,0,0,0,0,0,pbc*rr_bc_ooph,poc*rr_oc_ooph,0,0,psd,
        0,0,C,0,0,0,0,0,0,0,pbc*rr_bc_duet,poc*rr_oc_duet,0,0,psd,
        int_plp,int_plp,int_plp,C,0,0,0,0,0,0,pbc*rr1,poc*rr1,0,0,psd_cap_1,
        int_vus,int_vus,int_vus,0,C,0,0,0,0,0,pbc*rr2,poc*rr2,0,0,psd_cap_2,
        int_plp,int_plp,int_plp,0,0,C,0,0,0,0,pbc*rr4,poc*rr4,0,0,psd_cap_4,
        int_vus,int_vus,int_vus,0,0,0,C,0,0,0,pbc*rr5,poc*rr5,0,0,psd_cap_5,
        int_plp,int_plp,int_plp,0,0,0,0,C,0,0,pbc*rr7,poc*rr7,0,0,psd_cap_7,
        int_vus,int_vus,int_vus,0,0,0,0,0,C,0,pbc*rr8,poc*rr8,0,0,psd_cap_8,
        0,0,0,0,0,0,0,0,0,C,pbc,poc,0,0,psd,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pbcd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pocd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,0,ppostd,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        
        state_names = states
)
mat_none <- define_transition(
        C,0,0,0,0,0,0,0,0,0,pbc*rr_bc_mast,poc*rr_oc_mast,0,0,psd,
        0,C,0,0,0,0,0,0,0,0,pbc*rr_bc_ooph,poc*rr_oc_ooph,0,0,psd,
        0,0,C,0,0,0,0,0,0,0,pbc*rr_bc_duet,poc*rr_oc_duet,0,0,psd,
        0,0,0,C,0,0,0,0,0,0,pbc*rr1,poc*rr1,0,0,psd,
        0,0,0,0,C,0,0,0,0,0,pbc*rr2,poc*rr2,0,0,psd,
        0,0,0,0,0,C,0,0,0,0,pbc*rr4,poc*rr4,0,0,psd,
        0,0,0,0,0,0,C,0,0,0,pbc*rr5,poc*rr5,0,0,psd,
        0,0,0,0,0,0,0,C,0,0,pbc*rr7,poc*rr7,0,0,psd,
        0,0,0,0,0,0,0,0,C,0,pbc*rr8,poc*rr8,0,0,psd,
        0,0,0,0,0,0,0,0,0,C,pbc,poc,0,0,psd,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pbcd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,pocd,0,
        0,0,0,0,0,0,0,0,0,0,0,0,C,0,ppostd,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
        
        state_names = states
)

mat_none
plot(mat_none)

r <- 0.03
state_mast <- define_state(cost=0,QALY=discount(uint,r,first = TRUE))
state_ooph <- define_state(cost=0,QALY=discount(uint,r,first = TRUE))
state_duet <- define_state(cost=0,QALY=discount(uint,r,first = TRUE))
#this part was messed up in excel versione (revised version as comments)
# state_hp_plp <- define_state(cost=discount(0.5*(cmast+cooph),r),QALY=discount(unone-0.25*0.2,r))
# state_hp_vus <- define_state(cost=discount(0.1*(cmast+cooph),r),QALY=discount(unone-0.05*0.2,r))
# state_mp_plp <- define_state(cost=discount(0.5*(cmast+cooph),r),QALY=discount(unone-0.25*0.2,r))
# state_mp_vus <- define_state(cost=discount(0.1*(cmast+cooph),r),QALY=discount(unone-0.05*0.2,r))
# state_lp_plp <- define_state(cost=discount(0.5*(cmast+cooph),r),QALY=discount(unone-0.25*0.2,r))
# state_lp_vus <- define_state(cost=discount(0.1*(cmast+cooph),r),QALY=discount(unone-0.05*0.2,r))
state_hp_plp <- define_state(cost=discount(0.5*cmast,r,first = TRUE),QALY=discount(unone-0.25*0.2,r,first = TRUE))
state_hp_vus <- define_state(cost=discount(0.1,r,first = TRUE),QALY=discount(unone-0.05*0.2,r,first = TRUE))
state_mp_plp <- define_state(cost=discount(0.5*cmast,r,first = TRUE),QALY=discount(unone-0.25*0.2,r,first = TRUE))
state_mp_vus <- define_state(cost=discount(0.1,r,first = TRUE),QALY=discount(unone-0.05*0.2,r,first = TRUE))
state_lp_plp <- define_state(cost=discount(0.5*cmast,r,first = TRUE),QALY=discount(unone-0.25*0.2,r,first = TRUE))
state_lp_vus <- define_state(cost=discount(0.1,r,first = TRUE),QALY=discount(unone-0.05*0.2,r,first = TRUE))

state_np <- define_state(cost=0,QALY=discount(unone,r,first = TRUE))
state_breast <- define_state(cost=discount(cbc,r,first = TRUE),QALY=discount(ucancer,r,first = TRUE))
state_ovarian <- define_state(cost=discount(coc,r,first = TRUE),QALY=discount(ucancer,r,first = TRUE))
state_postcancer <- define_state(cost=discount(cpostc,r,first = TRUE),QALY=discount(upostc,r,first = TRUE))
state_cancerd <- define_state(cost=discount(cdeath,r,first = TRUE),QALY=0) #messed up in excel, not sure how to replicate here
# state_dead <- define_state(cost=discount(as.integer(state_time==1)*cdeath,r),QALY=0)
state_dead <- define_state(cost=0,QALY=0)

strat_none <- define_strategy(
        transition = mat_none,
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

strat_test <- define_strategy(
        transition = mat_test,
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
init_none <- c(mast=.01*1000L,
               ooph=.01*1000L,
               duet=.01*1000L,
               hp_plp=.97*1000L*gene_prev[1],
               hp_vus=.97*1000L*gene_prev[2],
               mp_plp=.97*1000L*gene_prev[4],
               mp_vus=.97*1000L*gene_prev[5],
               lp_plp=.97*1000L*gene_prev[7],
               lp_vus=.97*1000L*gene_prev[8],
               np=.97*1000L*sum(gene_prev[c(3,6,9,10)]))

init_test <- c(mast=1000L*(.05*sum(gene_prev[c(1,4,7)])+.03*sum(gene_prev[c(2,5,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
               ooph=1000L*(.03*sum(gene_prev[c(1,4,7)])+.02*sum(gene_prev[c(2,5,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
               duet=1000L*(.02*sum(gene_prev[c(1,2,4,5,7,8)])+.01*sum(gene_prev[c(3,6,9,10)])),
               hp_plp=.9*1000L*gene_prev[1],
               hp_vus=.93*1000L*gene_prev[2],
               mp_plp=.9*1000L*gene_prev[4],
               mp_vus=.93*1000L*gene_prev[5],
               lp_plp=.9*1000L*gene_prev[7],
               lp_vus=.93*1000L*gene_prev[8],
               np=.97*1000L*(gene_prev[3]+gene_prev[6]+gene_prev[9]+gene_prev[10]))

res_mod_none <- run_model(
        none = strat_none,
        parameters = param,
        cycles = 55,
        cost = cost,
        effect = QALY,
        state_time_limit = 1,
        init = c(as.numeric(init_none),0,0,0,0,0),
        method = "beginning"
)


res_mod_test <- run_model(
        test = strat_test,
        parameters = param,
        cycles = 55,
        cost = cost,
        effect = QALY,
        state_time_limit = 1,
        init = c(as.numeric(init_test),0,0,0,0,0),
        method = "beginning"
)

ct1 <- res_mod_none$eval_strategy_list$none$counts
addq1 <- 0.95*sum(init_none)-0.05*sum(init_none[c(4,6,8)])-0.01*sum(init_none[c(5,7,9)])
addc1 <- 12596*sum(init_none[c(4,6,8)])*0.5+sum(init_none[c(5,7,9)])*0.1
re1 <- res_mod_none$run_model %>% select(cost,QALY) %>% mutate(QALY=QALY+addq1,cost=cost+addc1)

ct2 <- res_mod_test$eval_strategy_list$test$counts
addq2 <- 0.95*sum(init_test)-0.05*sum(init_test[c(4,6,8)])-0.01*sum(init_test[c(5,7,9)])
addc2 <- 12596*sum(init_test[c(4,6,8)])*0.5+sum(init_test[c(5,7,9)])*0.1+658000
re2 <- res_mod_test$run_model %>% select(cost,QALY) %>% mutate(QALY=QALY+addq2, cost=cost+addc2)

comp <- rbind(re1,re2) %>% mutate(strategy=c("none","test"),ICER=(cost[2]-cost[1])/(QALY[2]-QALY[1])) #add test cost

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
