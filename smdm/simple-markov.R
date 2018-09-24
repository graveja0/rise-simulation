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
    disuA    = params$d_a/params$interval,
    costDrug = 365*params$c_tx/params$interval,
    costAlt  = 365*params$c_alt/params$interval,
    costTest = params$p_o * params$c_t,
    costBS   = params$c_bs,
    disuB    = params$d_b/params$interval,
    costBD   = params$c_bd,
    
    age_init = 40,
    t1 = model_time/params$interval, #model_time starts with 1
    t0 = (model_time-1)/params$interval,
    
    pD = map2_dbl(t0,t1,~gompertz_ratio(t0=.x,t1=.y,shape=params$shape,rate=params$rate)),
    pA = map_dbl(pD,~cap_max(value=rate_to_prob(params$r_a,1/params$interval),sd=.x)),
    pB = map_dbl(pD,~cap_max(value=rate_to_prob(params$r_b,1/params$interval),sd=.x)),
    
    fatalB = params$p_bd,
    pBS = pB*(1-fatalB),
    pBD = pB*fatalB,
    
    gene = 1, #0 or 1
    rr = map_dbl(gene, function(x) ifelse(x==0, 1, params$rr_b)),
    
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
  
  ####################
  # Healthy
  state_H <- define_state(
    cost     = 0,
    QALY     = discount(1/params$interval,dr),
    
    # Diagnostics
    A_acc    = 0,
    A_du_acc = 0,
    living   = 1,
    possible = discount(1/params$interval, dr),
    fatal_b  = 0,
    cost_g   = 0,
    cost_d   = 0,
    cost_tx  = 0,
    dis_a    = 0,
    dis_b    = 0
  )
  
  ####################
  # Indication 
  state_A <- define_state(
    cost = discount(dispatch_strategy(
      reference=ifelse(state_time==1,costA,0)+costDrug,
      genotype =ifelse(state_time==1,costA,0)+ifelse(state_time==1,costTest,0)+cDgenotype
    ), dr),
    QALY   = discount(ifelse(state_time<=params$d_at*params$interval,1/params$interval - disuA,1/params$interval),dr),
    
    # Diagnostics
    A_acc    = ifelse(state_time==1,1,0),
    A_du_acc = ifelse(state_time<=params$d_at*params$interval,1,0), # How many in disutility state
    living   = 1,
    possible = discount(1/params$interval, dr),
    fatal_b  = 0,
    cost_g   = discount(ifelse(state_time==1,costTest,0), dr),
    cost_d   = discount(dispatch_strategy(reference=costDrug,genotype=cDgenotype), dr),
    cost_tx  = discount(ifelse(state_time==1,costA,0), dr),
    dis_a    = discount(ifelse(state_time<=params$interval*params$d_at,disuA,0), dr),
    dis_b    = 0
  )
  
  ####################
  # Adverse Event Survivor
  state_BS <- define_state(
    cost = discount(dispatch_strategy(
      reference=ifelse(state_time==1,costBS,0)+costDrug,
      genotype=ifelse(state_time==1,costBS,0)+cDgenotype
    ), dr),
    QALY = discount(1/params$interval-disuB,dr),
    
    # Diagnostics
    A_acc    = 0,
    A_du_acc = 0,
    living   = 1,
    possible = discount(1/params$interval, dr),
    fatal_b  = 0,
    cost_g   = 0,
    cost_d   = discount(dispatch_strategy(reference=costDrug,genotype=cDgenotype), dr),
    cost_tx  = discount(ifelse(state_time==1,costBS,0), dr),
    dis_a    = 0,
    dis_b    = discount(disuB, dr)
  )

  ####################
  # Adverse Event Death
  state_BD <- define_state(
    cost = discount(ifelse(state_time==1,costBD,0),dr),
    QALY = 0,
    
    # Diagnostics
    A_acc    = 0,
    A_du_acc = 0,
    living   = 0,
    possible = 0,
    fatal_b  = ifelse(state_time==1,1,0),
    cost_g   = 0,
    cost_d   = 0,
    cost_tx  = discount(ifelse(state_time==1,costBD,0), dr),
    dis_a    = 0,
    dis_b    = 0
  )
  
  ####################
  # Secular Death
  state_D <- define_state(
    cost     = 0,
    QALY     = 0,
    
    # Diagnostics
    A_acc    = 0,
    A_du_acc = 0,
    living   = 0,
    possible = 0,
    fatal_b  = 0,
    cost_g   = 0,
    cost_d   = 0,
    cost_tx  = 0,
    dis_a    = 0,
    dis_b    = 0
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
    cycles = ceiling(params$horizon*params$interval),
    cost = cost,
    effect = QALY,
    state_time_limit=ceiling(min(params$interval*params$d_at+1, params$horizon*params$interval)),
    method="life-table"  # WTF? This should be called "trapezoidal" and better should be alternate simpsons extended!
  )
  
  ### add gene prevalence
  pop   <- data.frame(gene=c(0,1), .weights=c(100-params$p_g*100,params$p_g*100))
  res_h <- update(res_mod, newdata = pop)
  
  res_h
}

markov_summary <- function(solution, params)
{
  model <- if(params$p_o == 0.0) solution$combined_model$eval_strategy_list$reference else
           if(params$p_o == 1.0) solution$combined_model$eval_strategy_list$genotype  else
           stop("Probability of ordering test as a ratio not handled")

  summary <- if(params$p_o == 0.0) data.frame(solution$combined_model$run_model[1,]) else
                                   data.frame(solution$combined_model$run_model[2,])
  
  c(dCOST       = unname(summary[1,'cost']),
    dQALY       = unname(summary[1,'QALY']),
    possible    = unname(sum(model$values$possible)),
    fatal_b     = unname(sum(model$values$fatal_b)),
    living      = unname(model$values$living[length(model$values$living)]),
    disutil_a   = unname(sum(model$values$dis_a)),
    disutil_b   = unname(sum(model$values$dis_b)),
    dCOST.test  = unname(sum(model$values$cost_g)),
    dCOST.drug  = unname(sum(model$values$cost_d)),
    dCOST.treat = unname(sum(model$values$cost_tx))
  )/1000
}

markov_icer <- function(params, solution=NULL)
{
  if(is.null(solution)) solution   <- suppressMessages(markov_simulation(params))
  params$p_o <- 0
  reference  <- markov_summary(solution, params)
  params$p_o <- 1
  genotype   <- markov_summary(solution, params)

  c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
     NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
     dCOST.ref  = unname(reference['dCOST']),
     dCOST.test = unname(genotype['dCOST']),
     dQALY.ref  = unname(reference['dQALY']),
     dQALY.test = unname(genotype['dQALY'])
  )
}

#round(markov_summary(solution <- markov_simulation(params), params),5)
#markov_icer(params)