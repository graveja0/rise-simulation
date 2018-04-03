library(heemod)
library(diagram)
library(tidyverse)

rm(list=ls())

# parameters
load("raw_mortality.rda") #raw annual prob secular death, source data to fit death in simmer and numerical model.

doom <- function(value,cage) {
        ifelse(cage<max(lt$Age),value,0)
} # Once age reaches the maximum, pD becomes 1 and all other probs need to put 0.

doom2 <- function(value,cage) {
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
        age = age_init + markov_cycle,
        
        pD = map_dbl(age, function(x) doom2(lt$prob[lt$Age==x & lt$gender=="female"],x)),
        
        pA = map_dbl(age, function(x) doom(rescale_prob(p=0.1,from=10,to=1),x)),
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
        state_names = c("H","ABF","ABS","ABD","BS","BD","D"),
        C,pA*(1-pB),pA*pBS,pA*pBD,0,0,pD,
        0,C,0,0,pBS,pBD,pD,
        0,0,C,0,0,0,pD,
        0,0,0,1,0,0,0,
        0,0,0,0,C,0,pD,
        0,0,0,0,0,1,0,
        0,0,0,0,0,0,1
)

mat_genotype <- define_transition(
        state_names = c("H","ABF","ABS","ABD","BS","BD","D"),
        C,pA*(1-rr*pB),pA*pBS*rr,pA*pBD*rr,0,0,pD,
        0,C,0,0,rr*pBS,rr*pBD,pD,
        0,0,C,0,0,0,pD,
        0,0,0,1,0,0,0,
        0,0,0,0,C,0,pD,
        0,0,0,0,0,1,0,
        0,0,0,0,0,0,1
)


plot(mat_standard)

dr <- 0.03 # discounting rate
#states
state_H <- define_state(
        cost = 0,
        QALY = discount(1,dr)
)

state_ABF <- define_state(
        cost = discount(dispatch_strategy(
                standard=costA*ifelse(state_time<=1,1,0)+costDrug,
                genotype=costA*ifelse(state_time<=1,1,0)+cDgenotype
        ), dr),
        QALY = discount(1-disuA*ifelse(state_time<=1,1,0),dr)
)

state_ABS <- define_state(
        cost = discount(dispatch_strategy(
                standard=(costA+costBS)*ifelse(state_time<=1,1,0)+costDrug,
                genotype=(costA+costBS)*ifelse(state_time<=1,1,0)+cDgenotype
        ), dr),
        QALY = discount(1-disuB-disuA*ifelse(state_time<=1,1,0),dr)
)


state_BS <- define_state(
        cost = discount(dispatch_strategy(
                standard=costBS*ifelse(state_time<=1,1,0)+costDrug,
                genotype=costBS*ifelse(state_time<=1,1,0)+cDgenotype
        ), dr),
        QALY = discount(1-disuB,dr)
)

state_ABD <- define_state(
        cost = discount((costA+costBD)*ifelse(state_time<=1,1,0),dr),
        QALY = 0 
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
        ABF=state_ABF,
        ABS=state_ABS,
        ABD=state_ABD,
        BS=state_BS,
        BD=state_BD,
        D=state_D
)

strat_genotype <- define_strategy(
        transition = mat_genotype,
        H=state_H,
        ABF=state_ABF,
        ABS=state_ABS,
        ABD=state_ABD,
        BS=state_BS,
        BD=state_BD,
        D=state_D
)

# run
res_mod <- run_model(
        standard=strat_standard,
        genotype=strat_genotype,
        parameters = param,
        cycles = 70,
        cost = cost,
        effect = QALY,
        state_time_limit=1,
        method="life-table"
)


p <- data.frame(res_mod$eval_strategy_list$standard$parameters) #check parameters

res_mod$run_model
res_mod$eval_strategy_list$standard$counts
res_mod$eval_strategy_list$standard$values
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
res_h$combined_model$run_model %>% mutate(ICER=diff(cost)/diff(QALY))

plot(res_mod, type = "counts", panel = "by_state", free_y = TRUE,
     states=c("H","D","ABF","ABS","BS","ABD","BD")) +
        theme_bw() +
        scale_color_brewer(
                name = "Strategy",
                palette = "Set1"
        )

c1 <- res_h$combined_model$eval_strategy_list$standard$counts
c2 <- res_h$combined_model$eval_strategy_list$genotype$counts
v1 <- res_h$combined_model$eval_strategy_list$standard$values
v2 <- res_h$combined_model$eval_strategy_list$genotype$values


cp <- rbind(c1 %>% mutate(r=rownames(.),s=0), c2 %>% mutate(r=rownames(.),s=1))
cv <- rbind(v1 %>% mutate(r=rownames(.),s=0), v2 %>% mutate(r=rownames(.),s=1))

write.csv(cp,file="~/Desktop/M2yr_count.csv")
write.csv(cv,file="~/Desktop/M2yr_value.csv")
