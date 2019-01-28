library(heemod)
library(diagram)
library(tidyverse)
library(readxl)

###############################
# try to replicate BRCA model
###############################

### Comparison
# Results are almost identical
# Negative probability in the last cycle in no-testing scenario due competing intervention/death events

### Questions to ask/discuss
# computing costs for ovarian cancer deaths looks up the age in one cycle back
# post-both not at risk for breast cancer
# second surgery does not count in cost computation
# confusing QALY computation for interventions (I replicated this part outside of the simulation but was unable to interpret the formula, see function 'sp')

rm(list=ls()) 
###0. Read in params from excel
params <- read_excel(path = "~/Box Sync/*temp/RISE Model_BC_3.xlsm", sheet = "Breast_Ovarian", range = "B36:C91",col_types = c("text","numeric")) %>% set_names("param","value") %>%
        filter(!is.na(value))
pp <- unlist((params[,2]))
tb <- read_excel(path = "~/Box Sync/*temp/RISE Model_BC_3.xlsm", sheet = "Breast_Ovarian", range = "CX8:EL108",col_names=FALSE)

uptake_table <- tb %>% select(1,20,23) %>% set_names("age","mast","ooph")
uptake_table$none <- c(rep(as.numeric(params[14,2]),100),0)
bc_table <- tb %>% select(1,8,32,35) %>% set_names("age","no","mast","ooph")
oc_table <- tb %>% select(1,11,38,41) %>% set_names("age","no","mast","ooph")

life_table <- tb %>% select(1,2,4,5) %>% set_names("age","no","bc","oc")
#life_table[101,-1] <- life_table[101,-1]-0.000000000001

###1. Define inputs
#intervention
upt <- function(age,which) {
        uptake_table[uptake_table$age==age,][[which]]
}

#cancer
p_cancer <- function(age,cancer,which) {
        if(cancer=="bc") {
                p <- bc_table[bc_table$age==age,][[which]]
        }
        if(cancer=="oc") {
                p <- oc_table[oc_table$age==age,][[which]]
        }
        return(p)
}

#mortality
p_mort <- function(age,which) {
        life_table[life_table$age==age,][[which]]
}

param <- define_parameters(
        age_init = 45,
        age = age_init + model_time,
        
        #intervention uptakes depend on age
        upt_mast = map_dbl(age,~upt(age=.x,which="mast")),
        upt_ooph = map_dbl(age,~upt(age=.x,which="ooph")),
        upt_no = map_dbl(age,~upt(age=.x,which="none")),
        
        #cancer rates depend on age & intervention
        pbc_no = map_dbl(age,~p_cancer(cancer="bc",age=.x,which="no")),
        pbc_mast = map_dbl(age,~p_cancer(cancer="bc",age=.x,which="mast")),
        pbc_ooph = map_dbl(age,~p_cancer(cancer="bc",age=.x,which="ooph")),
        poc_no = map_dbl(age,~p_cancer(cancer="oc",age=.x,which="no")),
        poc_mast = map_dbl(age,~p_cancer(cancer="oc",age=.x,which="mast")),
        poc_ooph = map_dbl(age,~p_cancer(cancer="oc",age=.x,which="ooph")),
        
        #mortality depends on age & cancer status
        pd_no = map_dbl(age,~p_mort(age=.x,which="no")),
        pd_bc = map_dbl(age,~p_mort(age=.x,which="bc")),
        pd_oc = map_dbl(age,~p_mort(age=.x,which="oc")),
        
        uwell = 0.95,
        dknow = 0.05,
        
        dmast = 0.03,
        dooph = 0.03,
        
        ubc = 0.663,
        uoc = 0.628,
        upostbc = 0.81,
        upostoc = 0.72

)

###2. Define states & transition matrix (both scenarios have the same transition structure and initial states distribution, only differ in a few probs and testing cost)
states <- c("no_int","mast","ooph","post_mast","post_ooph","second_mast","second_ooph","post_both","bc1","bc2","oc1","oc2","dead","dead_bc","dead_oc")
mat_test <- define_transition(
        C,upt_mast,upt_ooph,0,0,0,0,0,pbc_no,0,poc_no,0,pd_no,0,0, #no_int
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, #mast
        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0, #ooph
        0,0,0,C,0,0,upt_ooph,0,pbc_mast,0,poc_mast,0,pd_no,0,0, #post_mast
        0,0,0,0,C,upt_mast,0,0,pbc_ooph,0,poc_ooph,0,pd_no,0,0, #post_ooph
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, #second_mast
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, #second_ooph
        0,0,0,0,0,0,0,C,0,0,poc_ooph,0,pd_no,0,0, #post_both
        0,0,0,0,0,0,0,0,0,C,0,0,0,pd_bc,0, #bc1
        0,0,0,0,0,0,0,0,0,C,0,0,0,pd_bc,0, #bc2
        0,0,0,0,0,0,0,0,0,0,0,C,0,0,pd_oc, #oc1
        0,0,0,0,0,0,0,0,0,0,0,C,0,0,pd_oc, #oc2
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead_bc
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead_oc
        
        state_names = states
)

mat_none <- define_transition(
        C,upt_no,upt_no,0,0,0,0,0,pbc_no,0,poc_no,0,pd_no,0,0, #no_int
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, #mast
        0,0,0,0,1,0,0,0,0,0,0,0,0,0,0, #ooph
        0,0,0,C,0,0,upt_ooph,0,pbc_mast,0,poc_mast,0,pd_no,0,0, #post_mast
        0,0,0,0,C,upt_mast,0,0,pbc_ooph,0,poc_ooph,0,pd_no,0,0, #post_ooph
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, #second_mast
        0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, #second_ooph
        0,0,0,0,0,0,0,C,0,0,poc_ooph,0,pd_no,0,0, #post_both
        0,0,0,0,0,0,0,0,0,C,0,0,0,pd_bc,0, #bc1
        0,0,0,0,0,0,0,0,0,C,0,0,0,pd_bc,0, #bc2
        0,0,0,0,0,0,0,0,0,0,0,C,0,0,pd_oc, #oc1
        0,0,0,0,0,0,0,0,0,0,0,C,0,0,pd_oc, #oc2
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead_bc
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0, #dead_oc
        
        state_names = states
)


mat_none
plot(mat_none)

###3. Define state cost/QALY implications
r <- 0.03
state_no_int <- define_state(cost=0,QALY=discount(uwell,r,first = TRUE))
#set QALY as 0, add back at the end
state_mast <- define_state(cost=discount(pp[24],r,first=TRUE),QALY=0)
state_ooph <- define_state(cost=discount(pp[25],r,first=TRUE),QALY=0) 
#state_mast <- define_state(cost=discount(pp[24],r,first=TRUE),QALY=discount(uwell-dmast,r,first = TRUE))
#state_ooph <- define_state(cost=discount(pp[25],r,first=TRUE),QALY=discount(uwell-dooph,r,first = TRUE))
state_post_mast <- define_state(cost=0,QALY=discount(uwell,r,first = TRUE))
state_post_ooph <- define_state(cost=0,QALY=discount(uwell,r,first = TRUE))
#no use
state_second_mast <- define_state(cost=0,QALY=0)
state_second_ooph <- define_state(cost=0,QALY=0)
# state_second_mast <- define_state(cost=discount(pp[24],r,first=TRUE),QALY=discount(uwell-dmast,r,first = TRUE))
# state_second_ooph <- define_state(cost=discount(pp[25],r,first=TRUE),QALY=discount(uwell-dooph,r,first = TRUE))
state_post_both <- define_state(cost=0,QALY=discount(uwell,r,first = TRUE))
state_bc1 <- define_state(cost=discount(pp[26]-as.integer(model_time>=20)*(pp[26]-pp[32]),r,first = TRUE),QALY=discount(ubc,r,first = TRUE))
state_bc2 <- define_state(cost=discount(pp[27],r,first = TRUE),QALY=discount(upostbc,r,first = TRUE))
state_oc1 <- define_state(cost=discount(pp[29]-as.integer(model_time>=20)*(pp[29]-pp[35]),r,first = TRUE),QALY=discount(uoc,r,first = TRUE))
state_oc2 <- define_state(cost=discount(pp[30],r,first = TRUE),QALY=discount(upostoc,r,first = TRUE))
state_dead <- define_state(cost=0,QALY=0)
state_dead_bc <- define_state(cost=discount(pp[28]-as.integer(model_time>=20)*(pp[28]-pp[34]),r,first = TRUE),QALY=0)
# typo here looking up based on the age in last cycle
state_dead_oc <- define_state(cost=discount(pp[31]-as.integer(model_time>=21)*(pp[31]-pp[37]),r,first = TRUE),QALY=0)

# sum(c(0,pp[c(24,25)],0,0,pp[c(24,25)],0,pp[c(26,27,29,30)],0,pp[c(28,31)])*ct1[1,])/1.03
# sum(c(0.95,0.92,0.92,0.95,0.95,0.92,0.92,0.95,0.663,0.81,0.628,0.72,0,0,0)*ct1[1,])/1.03

###4. Set up simulation
strat_none <- define_strategy(
        transition = mat_none,
        no_int=state_no_int,
        mast=state_mast,
        ooph=state_ooph,
        post_mast=state_post_mast,
        post_ooph=state_post_ooph,
        second_mast=state_second_mast,
        second_ooph=state_second_ooph,
        post_both=state_post_both,
        bc1=state_bc1,
        bc2=state_bc2,
        oc1=state_oc1,
        oc2=state_oc2,
        dead=state_dead,
        dead_bc=state_dead_bc,
        dead_oc=state_dead_oc
)

strat_test <- define_strategy(
        transition = mat_test,
        no_int=state_no_int,
        mast=state_mast,
        ooph=state_ooph,
        post_mast=state_post_mast,
        post_ooph=state_post_ooph,
        second_mast=state_second_mast,
        second_ooph=state_second_ooph,
        post_both=state_post_both,
        bc1=state_bc1,
        bc2=state_bc2,
        oc1=state_oc1,
        oc2=state_oc2,
        dead=state_dead,
        dead_bc=state_dead_bc,
        dead_oc=state_dead_oc
)

#load initial states
init_upt <- sum(tb[46,c(26,29)])
init_none <- define_init(no_int=1000L*(pp[3]*(1-pp[5])+pp[3]*pp[5]*(1-init_upt*3/2)),
                       mast=0,
                       ooph=0,
                       post_mast=1000L*(pp[3]*pp[5]*tb[46,26]),
                       post_ooph=1000L*(pp[3]*pp[5]*tb[46,29]),
                       second_mast=0,
                       second_ooph=0,
                       post_both=1000L*(pp[3]*pp[5]*init_upt/2),
                       bc1=0,
                       bc2=1000L*(tb[46,14]*pp[4]),
                       oc1=0,
                       oc2=1000L*(tb[46,17]*(1-pp[4])),
                       dead=0,
                       dead_bc=0,
                       dead_oc=0
)

###5. Run and summarize results
res_mod <- run_model(
        none = strat_none,
        test = strat_test,
        parameters = param,
        cycles = 55,
        cost = cost,
        effect = QALY,
        init = init_none,
        method = "beginning"
)

options(digits = 10)
ct1 <- res_mod$eval_strategy_list$none$counts
ct2 <- res_mod$eval_strategy_list$test$counts
addq1 <- sum(res_mod$eval_strategy_list$none$e_init*c(0.95,rep(0.92,2),rep(0.95,5),0.66,0.81,0.63,0.72,0,0,0))
addq2 <- sum(res_mod$eval_strategy_list$none$e_init*c(0.9,rep(0.92,2),rep(0.95,5),0.66,0.81,0.63,0.72,0,0,0))
addc <- sum(res_mod$eval_strategy_list$none$e_init*c(0,pp[24],pp[25],rep(0,2),pp[24],pp[25],0,pp[26],pp[27],pp[29],pp[30],0,pp[28],pp[31]))

#add back QALYs for intervention
sp <- function(dt) {
        dd <- dt$counts[,c("mast","ooph","post_ooph")]/1000
        dd$cyc <- seq(1:nrow(dd))
        dd$age <- seq(1:nrow(dd))+45
        dd$up_mast <- uptake_table$mast[uptake_table$age>45]
        dd$up_ooph <- uptake_table$ooph[uptake_table$age>45]
        dd <- dd %>% mutate(
                p1=mast*0.92,
                p2=ooph*(0.92+post_ooph*(0.95-up_mast*0.03-up_ooph*0.03)),
                #p3=0.92+post_ooph*(0.95-up_mast*0.03-up_ooph*0.03),
                dd=p1+p2,
                val=dd/(1+0.03)^cyc
        )
        1000*sum(dd$val)
}
addqq1 <- sp(dt=res_mod$eval_strategy_list$none)
addqq2 <- sp(dt=res_mod$eval_strategy_list$test)

res_mod$run_model %>% select(cost,QALY) %>% mutate(cost=(cost+c(addc,addc+657.5*1000))/1000,QALY=(QALY+c(addq1,addq2)+c(addqq1,addqq2))/1000,
                                                   strategy=c("None","Test"),ICER=c(NA,(cost[2]-cost[1])/(QALY[2]-QALY[1]))) 

# cc <- ct2
# cc$cyc <- seq(1:nrow(cc))
# p1 <- cc %>% mutate(vv=no_int*0.95/(1+0.03)^cyc) %>% pull(vv)/1000
# p2 <- cc %>% mutate(vv=(post_mast+post_ooph+post_both)*0.95/(1+0.03)^cyc) %>% pull(vv)/1000
# p3 <- cc %>% mutate(vv=bc1*0.663/(1+0.03)^cyc) %>% pull(vv)/1000
# p4 <- cc %>% mutate(vv=oc1*0.628/(1+0.03)^cyc) %>% pull(vv)/1000
# p5 <- cc %>% mutate(vv=(bc2*0.81+oc2*0.72)/(1+0.03)^cyc) %>% pull(vv)/1000
