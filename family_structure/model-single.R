library(deSolve)

if (!exists("ss_death")) ss_death <- read.csv("numerical/ss-death-2011.csv")

inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

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

params <- c(
  
  # Probabilities
  p_rt = 0.0,                  # Reactive Test
  p_pv = 0.0392,               # Pathogenic Variant Prevalence
  p_i1 = 0.05,                 # Probability of intervention one (Mastectomy)
  p_i2 = 0.03 ,                # Probability of intervention two (Oophrectomy)
  p_ib = 0.02,                 # Probability of both interventions
  
  # Risk Rates
  r_pt = 0, #inst_rate(0.2, 1), # Preemptive Test
  r_c  = 0.0012,               # Wild type, age dependent (see spreadsheet, SEER)
  r_dc = inst_rate(0.05,1),    # Risk of death from condition
  
  # Relative Risks
  rr_i1 = 0.1,                 # Intervention One relative risk
  rr_i2 = 0.51,                # Intervention Two relative risk
  rr_ib = 0.05,                # Intervention Both relative risk
  rr_pv = 15,                  # Pathegenic variant risk

  # Costs
  c_oc   = 75873,              # One time cost of condition
  c_dc   = 0,                  # Daily cost of condition
  c_oi1  = 12596,              # One time cost of intervention one
  c_oi2  =  8144,              # One time cost of intervention two
  c_oib  = 20740,              # One time cost of both interventions
  c_di1  = 0,                  # Daily cost of intervention one
  c_di2  = 0,                  # Daily cost of intervention two
  c_dib  = 0,                  # Daily cost of both interventions
  c_t    = 100,                # Cost of test

  # Disutility
  d_a1   = 0.15,               # Permanent Disutility for condition
  d_i1   = 0.05,               # Disutility for intervention one
  d_i2   = 0.05,               # Disutility for intervention two
  d_ib   = 0.10,               # Disutility for both interventions

  disc_rate =  inst_rate(1-1/1.03, 1)           # For computing discount
)

###################################
# Numerical Delay Differential Equation
genModel <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    
    r_d <- f_40yr_drate(t)

    list(c(
      uhw  = (-r_d -r_pt -r_c) * uhw,
      thw  = r_pt*uhw + (-r_d -r_c)*thw,
      ucw  = r_c*uhw + (-r_d -r_dc)*ucw,
      tcw  = r_c*thw+ (-r_d -r_dc)*tcw,
      uhp  = (-r_pt -r_d -r_c)*uhp,
      thp  = r_pt * (1-p_i1-p_i2-p_ib) * uhp + (-r_d-r_c) * thp,
      t1p  = r_pt * p_i1 * uhp + (-r_d-r_c*rr_i1) * t1p,
      t2p  = r_pt * p_i2 * uhp + (-r_d-r_c*rr_i2) * t2p,
      tbp  = r_pt * p_ib * uhp + (-r_d-r_c*rr_ib) * tbp,
      ucp  = r_c*uhp +(-r_d-r_dc)*ucp,
      tcp  = r_c*uhp + r_c*(thp+rr_i1*t1p+rr_i2*t2p+rr_ib*tbp)+ (-r_d-r_dc)*tcp,
      tacc = r_pt*(uhw+uhp),
      cdth = r_dc*(ucw+tcw+ucp+tcp),
      sdth = r_d*(uhw+thw+ucw+tcw+uhp+thp+t1p+t2p+tbp+ucp+tcp),
      disc = -disc_rate
    ))
  })
}

yinit <- with(as.list(c(params)), { c(
  uhw = 1-p_pv,
  thw = 0,
  ucw = 0,
  tcw = 0,
  uhp = p_pv,
  thp = 0,
  t1p = 0,
  t2p = 0,
  tbp = 0,
  ucp = 0,
  tcp = 0,
  tacc = 0,
  cdth = 0,
  sdth = 0,
  disc = 1
)})

# Important-- each time step is a day!!!
times <- seq(0, 10, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- dede(yinit, times, genModel, params)))

plot(out)

costs <- function(solution, params)
{
  k        <- length(solution[,1])
  simpson  <- alt_simp_coef(k)
  step     <- solution[2,'time'] - solution[1,'time']
  
  with(as.list(params), {

    # Daily costs work because timestep is a day
    test.cost <- c_t*(solution[1,'tacc']+sum(diff(solution[,'tacc'])*solution[2:k,'disc']))
    i1.cost   <- c_oi1*(solution[,'disc'] * solution[,'t1p'] )
    i2.cost   <- c_oi2*(solution[,'disc'] * solution[,'t2p'] )
    ib.cost   <- c_oib*(solution[,'disc'] * solution[,'tbp'] )
    intv.cost <- i1.cost + i2.cost + ib.cost
    ucw.cost  <- c_oc*(sum(diff(solution[,'ucw'])*solution[2:k,'disc'])) + 
                 c_dc*sum(solution[,'ucw']*solution[2:k,'disc'])
    tcw.cost  <- c_oc*(sum(diff(solution[,'tcw'])*solution[2:k,'disc'])) + 
                 c_dc*sum(solution[,'tcw']*solution[2:k,'disc'])
    ucp.cost  <- c_oc*(sum(diff(solution[,'ucp'])*solution[2:k,'disc'])) + 
                 c_dc*sum(solution[,'ucp']*solution[2:k,'disc'])
    tcp.cost  <- c_oc*(sum(diff(solution[,'tcp'])*solution[2:k,'disc'])) + 
                 c_dc*sum(solution[,'tcp']*solution[2:k,'disc'])
    cond.cost <- ucw.cost + tcw.cost + ucp.cost + tcp.cost
    
    c(dCOST       = unname(test.cost + intv.cost + cond.cost),
      
      # dQALY       = unname(pQALY-disA-disB),
      # possible    = unname(pQALY),
      # fatal_b     = unname(b_d),
      # living      = unname(life[k]),
      # disutil_a   = unname(disA),
      # disutil_b   = unname(disB),
      # 
      dCOST.test  = unname(test.cost),
      dCOST.intv  = unname(intv.cost),
      dCOST.cond  = unname(cond.cost)
     )
  })
}

# costs <- function(solution, params)
# {
#   k        <- length(solution[,1])
#   simpson  <- alt_simp_coef(k)
#   step     <- solution[2,'time'] - solution[1,'time']
#   
#   with(as.list(params), {
#     disc  <- length(key) * n + 3 # Discount Rate
#     tests <- disc - 1
# 
# 
#     # Testing costs
#     test.cost <- c_t*(solution[1,tests]+sum(diff(solution[,tests])*solution[2:k,disc]))
# 
#     # Loop over conditions
#     treatment.cost <- 0
#     drug.cost      <- 0
#     b_d            <- 0
#     for(i in 1:n)
#     {
#       # Compute Discounted Cost
#       treatment.cost <- treatment.cost +
#         c_a[i]  *sum(diff(solution[,maps('a_c',i)])*solution[2:k,disc]) +
#         c_bs[i] *sum(diff(solution[,maps('b_c',i)])*solution[2:k,disc]) +
#         c_bd[i] *sum(diff(solution[,maps('b_d',i)])*solution[2:k,disc])
#       drug.cost <- drug.cost + 
#         c_tx[i] *365*sum(simpson*solution[,maps('a_p',i)]*solution[,disc])*step +
#         c_alt[i]*365*sum(simpson*solution[,maps('a_a',i)]*solution[,disc])*step +
#         c_tx[i] *365*sum(simpson*solution[,maps('b_p',i)]*solution[,disc])*step +
#         c_alt[i]*365*sum(simpson*solution[,maps('b_a',i)]*solution[,disc])*step
#       
#       # Sum fatal b events from all conditions
#       b_d  <- b_d + solution[k,maps('b_d', i)]
#     }
#     
#     # Total living in model
#     life <- solution[,maps('h_u',1)] +
#             solution[,maps('h_t',1)] +
#             solution[,maps('a_p',1)] +
#             solution[,maps('a_a',1)] +
#             solution[,maps('b_p',1)] + 
#             solution[,maps('b_a',1)] 
# 
#     # Total possible life units is integral of discounted time
#     pQALY <- sum(simpson*life*solution[,disc])*step
# 
#     # Temp disutility of A
#     disA <- 0
#     disB <- 0
#     for(i in 1:n)
#     {
#       disA <- disA + d_a[i]*sum(simpson*solution[,maps('a_q',i)]*solution[,disc])*step
#       
#       # Permanent disutility for B (integration)
#       disB <- disB + 
#               d_b[i]*sum(simpson*solution[,maps('b_p',i)]*solution[,disc])*step + 
#               d_b[i]*sum(simpson*solution[,maps('b_a',i)]*solution[,disc])*step 
#     }
# 
#     c(dCOST       = unname(treatment.cost+test.cost+drug.cost),
#       dQALY       = unname(pQALY-disA-disB),
#       possible    = unname(pQALY),
#       fatal_b     = unname(b_d),
#       living      = unname(life[k]),
#       disutil_a   = unname(disA),
#       disutil_b   = unname(disB),
#       dCOST.test  = unname(test.cost),
#       dCOST.drug  = unname(drug.cost),
#       dCOST.treat = unname(treatment.cost)
#       )
#   })
# }
# 
# # Defined scenarios
# scenarios <- c("none", "reactive-single", "reactive-panel", "preemptive-panel")
# 
# generate.params <- function(config, i, scenario, disc_rate = inst_rate(0.03, 1))
# {
#   risks        <- unlist(config$risk[i,])
#   disutilities <- unlist(config$disutility[i,])
#   durations    <- unlist(config$duration[i,])
#   names(durations) <- names(config$duration)
#   costs        <- unlist(config$cost[i,])
#   
#   # Start building the params list with the length
#   n            <- length(unique(gsub(".*_SC_","",names(risks[grep("_SC_",names(risks))]))))
#   params       <- list(n=n)
#   
#   if(n != round(n)) stop("Error (n) In generate.params")
#   
#   getval <- function(x,tt) unname(tt[grep(x,names(tt))])
#   
#   params$p_p   <- if(scenario == "reactive-panel") 1.0 else 0.0
#   params$p_o   <- if(scenario == "none") rep(0.0, n) else getval("vProbabilityOrder",risks)
#   params$p_r   <- if(scenario == "none") rep(0.0, n) else getval("vProbabilityRead",risks)
#   params$p_bd  <- getval("vFatalB_",risks) # Probability of death as direct result of B
#   params$p_g   <- getval("vGene_",risks) # Probability of genetic variant
#   params$r_a   <- unname(inst_rate(getval("vRiskA_",risks), getval("vDurationA_",risks))) # Rate of a
#   params$r_b   <- unname(inst_rate(getval("vRiskB_",risks), getval("vDurationB_",risks))) # Rate of b
#   params$rr_b  <- getval("vRR_B_",risks) # Relative Risk of B when on alt treatment
# 
#   # Costs
#   params$c_a   <- getval("A_c_",costs)       # Cost of Event A
#   params$c_bs  <- getval("B_Survive_",costs) # Cost of Surviving Event B
#   params$c_bd  <- getval("B_Death_",costs )  # Cost of Death from Event B
#   params$c_tx  <- getval("rx_",costs)        # Cost of Treatment (Daily)
#   params$c_alt <- getval("alt_",costs)       # Cost of alternate treatment (Daily)
#   
#   params$c_t   <- if(scenario %in% c("reactive-panel", "preemptive-panel"))
#                   { config$global$panel_test[1] } else { config$global$single_test[1] }
#                   
#   params$d_a   <- getval("A_",disutilities)         # Disutility of A
#   params$d_at  <- getval("A_",durations)/365        # Duration of A in years.
#   params$d_b   <- getval("B_Survive_",disutilities) # Disutility of B
#   
#   params$disc_rate <- disc_rate       # For computing discount
# 
#   params
# }
# 
# generate.initial <- function(scenario, params)
# {
#   n <- params[["n"]]
#   yinit <- rep(0, length(key)*n+1)
#   yinit <- c(yinit, 1)
#   
#   if(scenario == "preemptive-panel")
#   {
#     for(i in 1:n) yinit[map("h_t", i)] <- 1
#     yinit[length(yinit)-1] <- 1 # Second to last is number of tests.
#   } else {
#     for(i in 1:n) yinit[map("h_u", i)] <- 1
#   }
#   
#   yinit
# }
# 
# # Config is latin hyper cube variable: drawn.parameter.values
# # i is the point to use from the cube
# # scenario specifies scenario
# # times is the time points to solve (resolution)
# model.run <- function(config, i, scenario, times=seq(0, 80, by=2/365))
# {
#   params  <- generate.params(config, i, scenario)
#   init    <- generate.initial(scenario, params)
#   x       <- dede(init, times, genModel, params)
#   costs(x, params)
# }

