library(deSolve)

# Vector must be 8 points or longer for alternate simpson integration
alt_ext_simpson <- function(x)
  sum(c(17/48, 59/48, 43/48, 49/48, rep(1, length(x)-8), 49/48, 43/48, 59/48, 17/48)*x)


###################################
# Numerical Delay Differential Equation
genModel <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    r_d <- rate*exp(shape*t) # Gompertz model of death from params
    
    # exp(- integral of death risk from t-d_at to t)
    death  <- if(t < d_at) 0 else exp(-rate*(exp(shape*t) - exp(shape*(t-d_at)))/shape)
    left_p <- if(t < d_at) 0 else lagderiv(t-d_at, 7)*exp(-r_b*d_at)*death
    left_a <- if(t < d_at) 0 else lagderiv(t-d_at, 8)*exp(-r_b*rr_b*d_at)*death
    
    list(c(
      h_u = (-r_a-r_d)*h_u,                            # Healthy Untested
      a_p = r_a*(1-p_o*p_g)*h_u-r_b*a_p  -r_d*a_p,     # Indication Primary Treat
      a_a = r_a*p_o*p_g*h_u-r_b*rr_b*a_a -r_d*a_a,     # Indication Alternate Treat
      a_c = r_a*h_u,                                   # Cumulative Indications
      a_dp = r_a*(1-p_o*p_g)*h_u-r_b*a_dp-r_d*a_dp-left_p, # Temporary Disutility of Indication Primary Occupation
      a_da = r_a*p_o*p_g*h_u-r_b*rr_b*a_da-r_d*a_da-left_a, # Temporary Disutility of Indication Alternate Occupation
      a_ep = r_a*(1-p_o*p_g)*h_u,                      # Helps compute leftover to clear from a_dp
      a_ea = r_a*p_o*p_g*h_u,                          # Helps compute leftover to clear from a_da
      b_p = r_b*(1-p_bd)*a_p-r_d*b_p,                  # Adverse Event Primary Treat
      b_a = r_b*rr_b*(1-p_bd)*a_a -r_d*b_a,            # Adverse Event Alternate Treat
      b_d = r_b*p_bd*a_p +r_b*rr_b*p_bd*a_a,           # Deaths from Adverse Event
      b_c = r_b*(1-p_bd)*a_p +r_b*rr_b*(1-p_bd)*a_a,   # Cumulative Adverse Events
      tests = p_o*r_a*h_u,                             # Cumulative Tests
      dsc = -disc_rate*dsc                             # Discount Curve
    ))
  })
}

deq_summary <- function(solution, params)
{
  k        <- length(solution[,1])
  step     <- solution[2,'time'] - solution[1,'time']
  
  with(as.list(params), {
    # Testing costs
    test.cost <- c_t*(solution[1,"tests"]+sum(diff(solution[,"tests"])*solution[2:k,"dsc"]))
    
    # Compute dscounted Cost
    treatment.cost <- 
      c_a  *alt_ext_simpson(diff(solution[,"a_c"])*solution[2:k,"dsc"]) +
      c_bs *alt_ext_simpson(diff(solution[,"b_c"])*solution[2:k,"dsc"]) +
      c_bd *alt_ext_simpson(diff(solution[,"b_d"])*solution[2:k,"dsc"])
    drug.cost <-
      c_tx *365*alt_ext_simpson(solution[,"a_p"]*solution[,"dsc"])*step +
      c_alt*365*alt_ext_simpson(solution[,"a_a"]*solution[,"dsc"])*step +
      c_tx *365*alt_ext_simpson(solution[,"b_p"]*solution[,"dsc"])*step +
      c_alt*365*alt_ext_simpson(solution[,"b_a"]*solution[,"dsc"])*step
    
    # Total living in model
    life <- rowSums(solution[,c("h_u","a_p","a_a","b_p","b_a")])
    
    # Total possible life units is integral of dscounted time
    pQALY <- alt_ext_simpson(life*solution[,"dsc"])*step
    
    # Temp Disutility of Indication
    disA <- d_a*alt_ext_simpson(solution[,'a_dp']*solution[,"dsc"])*step +
            d_a*alt_ext_simpson(solution[,'a_da']*solution[,"dsc"])*step
    # Permanent Disutility for Adverse Event
    disB <- d_b*alt_ext_simpson(solution[,'b_p']*solution[,"dsc"])*step + 
            d_b*alt_ext_simpson(solution[,'b_a']*solution[,"dsc"])*step 
    
    c(dCOST       = unname(treatment.cost+test.cost+drug.cost),
      dQALY       = unname(pQALY-disA-disB),
      possible    = unname(pQALY),
      fatal_b     = unname(solution[k,"b_d"]),
      living      = unname(life[k]),
      disutil_a   = unname(disA),
      disutil_b   = unname(disB),
      dCOST.test  = unname(test.cost),
      dCOST.drug  = unname(drug.cost),
      dCOST.treat = unname(treatment.cost)
    )
  })
}

deq_simulation <- function(params)
{
  params$disc_rate <- inst_rate(1-1/(1 + params$disc), 1)
  
  init  <- c(h_u=1, a_p=0, a_a=0, a_c=0, a_dp=0, a_da=0, a_ep=0, a_ea=0, b_p=0, b_a=0, b_d=0, b_c=0, tests=0, dsc=1)
  
  times <- seq(0, params$horizon, by=params$resolution)
  
  dede(init, times, genModel, params)
}

deq_icer <- function(params, reference=NULL, genotype=NULL)
{
  params$p_o <- 0.0 # No testing, reference
  reference <- deq_summary(if(is.null(reference)) deq_simulation(params) else reference, params)
           
  
  params$p_o <- 1.0 # Genotype testing upon indication
  genotype   <- deq_summary(if(is.null(genotype)) deq_simulation(params) else genotype, params)
  
  c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
     NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
     dCOST.ref  = unname(reference['dCOST']),
     dCOST.test = unname(genotype['dCOST']),
     dQALY.ref  = unname(reference['dQALY']),
     dQALY.test = unname(genotype['dQALY'])
  )
}
