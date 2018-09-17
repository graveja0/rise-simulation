library(heemod)
library(dplyr)
library(purrr)
library(flexsurv)

source("simple-params.R")

# secular death rate over time interval
gompertz_ratio <- function(t0, t1, shape, rate)
{
  r <- (pgompertz(t1, shape, rate) - pgompertz(t0, shape, rate)) / (1 - pgompertz(t0, shape, rate))
  if(is.na(r)) r=1
  return(r)
}

# Once secular death mortablity reaches 1, all other probs need to put 0. 
cap_max <- function(value,sd) ifelse(value+sd>1,1-sd,value)


markov_simulation <- function(params)
{
  param <- define_parameters(
    costA    = params$c_a,
    disuA    = params$d_a,
    costDrug = 365*params$c_tx/params$interval,
    costAlt  = 365*params$c_alt/params$interval,
    costBS   = params$c_bs,
    disuB    = params$d_b,
    costBD   = params$c_bd,
    
    age_init = 40,
    t1 = model_time/params$interval, #model_time starts with 1
    t0 = (model_time-1)/params$interval,
    
    pD = map2_dbl(t0,t1,~gompertz_ratio(t0=.x,t1=.y,shape=params$shape,rate=params$rate)),
    pA = map_dbl(pD,~cap_max(value=rescale_prob(p=0.1,from=10,to=1/params$interval),sd=.x)),
    pB = map_dbl(pD,~cap_max(value=rescale_prob(p=0.02,from=1,to=1/params$interval),sd=.x)),
    
    fatalB = params$p_bd,
    pBS = pB*(1-fatalB),
    pBD = pB*fatalB,
    
    gene = 1, #0 or 1
    rr = 1-params$rr_b*gene,
    
    #just for genotype strategy
    cDgenotype = map_dbl(gene, function(x) ifelse(x==0,costDrug,costAlt))
  )
  
  #transition matrix
  mat_reference <- define_transition(
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
  
  dr <- rescale_discount_rate(params$disc,from=params$interval,to=1) # discounting rate
  #states
  state_H <- define_state(
    cost = 0,
    QALY = discount(1,dr)
  )
  
  state_A <- define_state(
    cost = discount(dispatch_strategy(
      reference=costA*ifelse(state_time<=1,1,0)+costDrug,
      genotype=costA*ifelse(state_time<=1,1,0)+cDgenotype
    ), dr),
    QALY = discount(1-disuA*ifelse(state_time<=params$interval,1,0),dr)
  )
  
  
  state_BS <- define_state(
    cost = discount(dispatch_strategy(
      reference=costBS*ifelse(state_time<=1,1,0)+costDrug,
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
  strat_reference <- define_strategy(
    transition = mat_reference,
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
    reference=strat_reference,
    genotype=strat_genotype,
    parameters = param,
    cycles = params$horizon*params$interval,
    cost = cost,
    effect = QALY,
    state_time_limit=params$interval,
    method="life-table"  # WTF? This crap again? This should be alternate simpsons extended!
  )
  
  ### add gene prevalence
  pop   <- data.frame(gene=c(0,1), .weights=c(100-params$p_g*100,params$p_g*100))
  res_h <- update(res_mod, newdata = pop)
  
  # compare
  res_h
}

# Goal: dCOST    dQALY possible     fatal_b    living   disutil_a disutil_b dCOST.test dCOST.drug dCOST.treat
markov_summary <- function(solution, params)
{
  solution$combined_model$run_model %>%
  dplyr::mutate(ICER=diff(cost)/diff(QALY)*params$interval) %>%
  data.frame()
}

markov_icer <- function(params)
{
  ot <- markov_summary(markov_simulation(params), params)
  ot <- select(ot, strategy=.strategy_names,cost,QALY) %>%
        mutate(cost=cost/1000,QALY=QALY/(1000*params$interval))
  c(ICER=(ot$cost[2]-ot$cost[1])/(ot$QALY[2]-ot$QALY[1]),
    NMB=unname((ot$cost[1] - ot$cost[2]) + params$wtp*(ot$QALY[1] - ot$QALY[2])),
    dCOST.ref=ot$cost[1],
    dCOST.test=ot$cost[2],
    dQALY.ref=ot$QALY[1],
    dQALY.test=ot$QALY[2]
    )
}

markov_icer(params)