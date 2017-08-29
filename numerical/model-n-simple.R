library(deSolve)

if (!exists("ss_death")) ss_death <- read.csv("ss-death-2011.csv")

inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

#source('params-3.R')

###################################
# Numerical approach to secular death (very high accuracy!)
f_40yr_percent_d    <- c(ss_death$f_death_prob[41:120])
sim_adj_age         <- 0:79 + 0.5 # 0.5 offset since percentage is for whole year
f_40yr_per_d_spline <- splinefun(sim_adj_age, f_40yr_percent_d)
plot(1:2); dev.off()

# Clamped at infinite rate via pmin
f_40yr_drate <- function(t) inst_rate(pmin(f_40yr_per_d_spline(t), 1),1)

# This is for doing numberical integration of a set of numbers at an even interval
alt_simp_coef <- function(i)
{
  if (i < 8) stop("Invalid Simpson coefficient size")
  
  # pg.117 4.1.14, Numerical Recipes in C, 1st edition
  c(17/48, 50/48, 43/48, 49/48, rep(1, i-8), 49/48, 43/48, 50/48, 17/48) 
}

key <- list(
  h_u  =  1,
  h_t  =  2,
  a_p  =  3,
  a_a  =  4,
  a_c  =  5,
  a_l  =  6,
  a_e  =  7,
  a_q  =  8,
  b_p  =  9,
  b_a  = 10,
  b_d  = 11,
  b_c  = 12
)

map <- function(name, n) key[[name]] + (n-1) * length(key)
maps <- function(name, n) key[[name]] + (n-1) * length(key) + 1  # Solution is offset by time dimension

###################################
# Numerical Delay Differential Equation
genModel <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    
    # Interaction is primarily through death (other is panel test)
    # Living (should be same in all models i)
    liv  <- y[map("h_u", 1)] + 
            y[map("h_t", 1)] + 
            y[map("a_p", 1)] + 
            y[map("a_a", 1)] + 
            y[map("b_p", 1)] +
            y[map("b_a", 1)]
    
    rate <- c()
    for(i in 1:n)
    {
      r_d <- f_40yr_drate(t)
      for(j in 1:n) if(j != i) r_d <- r_d + (r_b[j]*p_bd[j]*(y[map("a_p", j)]+rr_b[j]*y[map("a_a", j)]))/liv
      
      r_p <- 0
      if(p_p > 0) for(j in 1:n) if(j != i) r_p <- r_p + p_p*p_o[i]*r_a[j]*y[map("h_u", j)] / liv

      rate <- c(rate, 
        (-r_p-r_a[i]-r_d)*y[map("h_u", i)],
        r_p*y[map("h_u", i)] + (-r_a[i]-r_d)*y[map("h_t", i)],
        r_a[i]*(1-p_o[i]*p_g[i])*y[map("h_u", i)]+r_a[i]*(1-p_g[i]*p_r[i])*y[map("h_t", i)]-r_b[i]*y[map("a_p", i)] -r_d*y[map("a_p", i)],
        r_a[i]*p_o[i]*p_g[i]*y[map("h_u", i)]+r_a[i]*p_g[i]*p_r[i]*y[map("h_t", i)] -r_b[i]*rr_b[i]*y[map("a_a", i)] -r_d*y[map("a_a", i)],
        r_a[i]*y[map("h_u", i)]+r_a[i]*y[map("h_t", i)],
        r_d*y[map("a_q", i)],
        r_d*y[map("a_q", i)] - if(t < d_at[i]) 0 else lagderiv(t-d_at[i], map("a_l",i)),
        r_a[i]*y[map("h_u", i)]+r_a[i]*y[map("h_t", i)] - if(t < d_at[i]) 0 else lagderiv(t-d_at[i], map("a_c", i))*exp(-y[map("a_e", 1)]),
        r_b[i]*(1-p_bd[i])*y[map("a_p", i)]-r_d*y[map("b_p", i)],
        r_b[i]*rr_b[i]*(1-p_bd[i])*y[map("a_a", i)] -r_d*y[map("b_a", i)],
        r_b[i]*p_bd[i]*y[map("a_p", i)] +r_b[i]*rr_b[i]*p_bd[i]*y[map("a_a", i)],
        r_b[i]*(1-p_bd[i])*y[map("a_p", i)] +r_b[i]*rr_b[i]*(1-p_bd[i])*y[map("a_a", i)]
      )
    }

    tests <- p_o[i]*sum(sapply(1:n, function(i) {
      r_a[i]*y[map("h_u", i)]
    }))
      
    list(c(rate,
           tests,
           -disc_rate*y[n * length(key) + 2]
    ))
  })
}


# 
# 
# times <- seq(0, 40, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
# print(system.time(out <- dede(yinit, times, genModel, params)))
# 
# plot(out)

costs <- function(solution, params)
{
  k        <- length(solution[,1])
  simpson  <- alt_simp_coef(k)
  step     <- solution[2,'time'] - solution[1,'time']
  
  with(as.list(params), {
    disc  <- length(key) * n + 3 # Discount Rate
    tests <- disc - 1


    # Testing costs
    test.cost <- c_t*(solution[1,tests]+sum(diff(solution[,tests])*solution[2:k,disc]))

    # Loop over conditions
    treatment.cost <- 0
    drug.cost      <- 0
    b_d            <- 0
    for(i in 1:n)
    {
      # Compute Discounted Cost
      treatment.cost <- 
        c_a[i]  *sum(diff(solution[,maps('a_c',i)])*solution[2:k,disc]) +
        c_bs[i] *sum(diff(solution[,maps('b_c',i)])*solution[2:k,disc]) +
        c_bd[i] *sum(diff(solution[,maps('b_d',i)])*solution[2:k,disc])
      drug.cost <-
        c_tx[i] *365*sum(simpson*solution[,maps('a_p',i)]*solution[,disc])*step +
        c_alt[i]*365*sum(simpson*solution[,maps('a_a',i)]*solution[,disc])*step +
        c_tx[i] *365*sum(simpson*solution[,maps('b_p',i)]*solution[,disc])*step +
        c_alt[i]*365*sum(simpson*solution[,maps('b_a',i)]*solution[,disc])*step
      
      # Sum fatal b events from all conditions
      b_d  <- b_d + solution[k,maps('b_d', i)]
    }
    
    # Total living in model
    life <- solution[,maps('h_u',1)] +
            solution[,maps('h_t',1)] +
            solution[,maps('a_p',1)] +
            solution[,maps('a_a',1)] +
            solution[,maps('b_p',1)] + 
            solution[,maps('b_a',1)] 

    # Total possible life units is integral of discounted time
    pQALY <- sum(simpson*life*solution[,disc])*step

    # Temp disutility of A
    disA <- 0
    disB <- 0
    for(i in 1:n)
    {
      disA <- disA + d_a[i]*sum(simpson*solution[,maps('a_q',i)]*solution[,disc])*step
      
      # Permanent disutility for B (integration)
      disB <- disB + 
              d_b[i]*sum(simpson*solution[,maps('b_p',i)]*solution[,disc])*step + 
              d_b[i]*sum(simpson*solution[,maps('b_a',i)]*solution[,disc])*step 
    }

    c(dCOST       = unname(treatment.cost+test.cost+drug.cost),
      dQALY       = unname(pQALY-disA-disB),
      possible    = unname(pQALY),
     # disutil_a  = unname(disA),
     # disutil_b  = unname(disB),
      fatal_b     = unname(b_d),
      living      = unname(life[k]),
      dCOST.test  = unname(test.cost),
      dCOST.drug  = unname(drug.cost),
      dCOST.treat = unname(treatment.cost)
      )
  })
}

# Defined scenarios
scenarios <- c("none", "reactive-single", "reactive-panel", "preemptive-panel")

generate.params <- function(config, i, scenario, disc_rate = inst_rate(0.03, 1))
{
  risks        <- unlist(config$risk[i,])
  disutilities <- unlist(config$disutility[i,])
  durations    <- unlist(config$duration[i,])
  names(durations) <- names(config$duration)
  costs        <- unlist(config$cost[i,])
  
  # Start building the params list with the length
  n            <- length(unique(gsub(".*_SC_","",names(risks[grep("_SC_",names(risks))]))))
  params       <- list(n=n)
  
  if(n != round(n)) stop("Error (n) In generate.params")
  
  getval <- function(x,tt) unname(tt[grep(x,names(tt))])
  
  params$p_p   <- if(scenario == "reactive-panel") 1.0 else 0.0
  params$p_o   <- if(scenario == "none") rep(0.0, n) else getval("vProbabilityOrder",risks)
  params$p_r   <- if(scenario == "none") rep(0.0, n) else getval("vProbabilityRead",risks)
  params$p_bd  <- getval("vFatalB_",risks) # Probability of death as direct result of B
  params$p_g   <- getval("vGene_",risks) # Probability of genetic variant
  params$r_a   <- unname(inst_rate(getval("vRiskA_",risks), getval("vDurationA_",risks))) # Rate of a
  params$r_b   <- unname(inst_rate(getval("vRiskB_",risks), getval("vDurationB_",risks))) # Rate of b
  params$rr_b  <- getval("vRR_B_",risks) # Relative Risk of B when on alt treatment

  # Costs
  params$c_a   <- getval("A_c_",costs)       # Cost of Event A
  params$c_bs  <- getval("B_Survive_",costs) # Cost of Surviving Event B
  params$c_bd  <- getval("B_Death_",costs )  # Cost of Death from Event B
  params$c_tx  <- getval("rx_",costs)   # Cost of Treatment (Daily)
  params$c_alt <- getval("alt_",costs)    # Cost of alternate treatment (Daily)
  
  params$c_t   <- if(scenario %in% c("reactive-panel", "preemptive-panel"))
                  { config$global$panel_test[1] } else { config$global$single_test[1] }
                  
  params$d_a   <- getval("A_",disutilities) # Disutility of A
  params$d_at  <- getval("A_",durations)/365    # Duration of A in years.
  params$d_b   <- getval("B_Survive_",disutilities) # Disutility of B
  
  params$disc_rate <- disc_rate       # For computing discount

  params
}

generate.initial <- function(scenario, params)
{
  n <- params[["n"]]
  yinit <- rep(0, length(key)*n+1)
  yinit <- c(yinit, 1)
  
  if(scenario == "preemptive-panel")
  {
    for(i in 1:n) yinit[map("h_t", i)] <- 1
    yinit[length(yinit)-1] <- 1 # Second to last is number of tests.
  } else {
    for(i in 1:n) yinit[map("h_u", i)] <- 1
  }
  
  yinit
}

# Config is latin hyper cube variable: drawn.parameter.values
# i is the point to use from the cube
# scenario specifies scenario
# times is the time points to solve (resolution)
model.run <- function(config, i, scenario, times=seq(0, 80, by=2/365))
{
  params  <- generate.params(config, i, scenario)
  init    <- generate.initial(scenario, params)
  x       <- dede(init, times, genModel, params)
  costs(x, params)
}

