library(heemod)
library(tidyverse)
library(flexsurv)

# parameters
load("raw_mortality.rda") #raw annual prob secular death, source data to fit death in simmer and numerical model.

#secular death rate
gompertz_ratio <- function(t0, t1, shape, rate)
{
        r <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) / (1 - pgompertz(t0, shape, rate))
        if(is.na(r)) r=1
        return(r)
}

# Gompertz model for 40yr female
shape <- 0.1007511
rate  <- 0.0008370717

# Once secular death mortablity reaches 1, all other probs need to put 0. 
cap_max <- function(value,sd) {
        ifelse(sd==1,0,value)
}

param <- define_parameters(
        costA = 10000,
        disuA = 0.05,
        costDrug = 365*0.5/interval,
        costAlt = 365*5/interval,
        costBS = 25000,
        disuB = 0.02,
        costBD = 15000,
        
        age_init = 40,
        t1 = model_time/interval, #model_time starts with 1
        t0 = (model_time-1)/interval,
        
        pD = map2_dbl(t0,t1,~gompertz_ratio(t0=.x,t1=.y,shape=shape,rate=rate)),
        pA = map_dbl(pD,~cap_max(value=rescale_prob(p=0.1,from=10,to=1/interval),sd=.x)),
        pB = map_dbl(pD,~cap_max(value=rescale_prob(p=0.02,from=1,to=1/interval),sd=.x)),
        
        fatalB = 0.05,
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

dr <- rescale_discount_rate(0.03,from=interval,to=1) # discounting rate
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
        QALY = discount(1-disuA*ifelse(state_time<=interval,1,0),dr)
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
        cycles = 85*interval,
        cost = cost,
        effect = QALY,
        state_time_limit=interval,
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
res_h$combined_model$run_model %>% mutate(ICER=diff(cost)/diff(QALY)*interval)

c1 <- res_h$combined_model$eval_strategy_list$standard$counts
c2 <- res_h$combined_model$eval_strategy_list$genotype$counts
v1 <- res_h$combined_model$eval_strategy_list$standard$values
v2 <- res_h$combined_model$eval_strategy_list$genotype$values

cp <- rbind(c1 %>% mutate(r=rownames(.),s=0), c2 %>% mutate(r=rownames(.),s=1))
cv <- rbind(v1 %>% mutate(r=rownames(.),s=0), v2 %>% mutate(r=rownames(.),s=1))
