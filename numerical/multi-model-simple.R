library(deSolve)
library(flexsurv) # For pgompertz

ss_death <- read.csv("ss-death-2011.csv")

inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

###################################
# Main model parameters
params <- c(

  # Probabilitieis
  p_bd1 = 0.05,               # Probability of death as direct result of B
  p_g1  = 0.2,                # Prevalence of targeted gene variant
  p_bd2 = 0.05,               # Probability of death as direct result of B
  p_g2  = 0.2,                # Prevalence of targeted gene variant

  p_o   = 0.75,               # Probability of ordering test
  p_r   = 0.75,               # Probability of reading test if it exists
  p_p   = 1,                  # Probability of panel test when test ordered

  # Risks
  r_a1   = inst_rate(0.7, 10), # Rate of healthy having event A over 10 years
  r_b1   = inst_rate(0.2, 3),  # Rate of post-A  having event B over 3 years
  rr_b1  = 0.5,                # Relative Risk for Event B with better treatment informed by genotyping
  r_a2   = inst_rate(0.2, 10), # Rate of healthy having event A over 10 years
  r_b2   = inst_rate(0.5, 3),  # Rate of post-A  having event B over 3 years
  rr_b2  = 0.5,                # Relative Risk for Event B with better treatment informed by genotyping
  
  # Costs
  c_a   = 10000,              # Cost of event A
  c_bs  = 25000,              # Cost for event B survival
  c_bd  = 15000,              # Cost for event B death
  c_tx  = 100/30,             # Daily cost of standard treatment
  c_alt = 300/30,             # Daily cots of alternate treatment
  c_t   = 100,                # Cost of test

  # Disutility
  d_a   = 0.2,                # Disutility for event A
  d_at  = 3*365,              # Time in days for disutility of event A
  d_b   = 0.15,               # Permanent Disutility for event B
  

  disc_rate = 1e-12           # For computing discount
)

###################################
# Numerical approach to secular death (very high accuracy!)
f_40yr_percent_d    <- c(ss_death$f_death_prob[41:120])
sim_adj_age         <- 0:79 + 0.5 # 0.5 offset since percentage is for whole year
f_40yr_per_d_spline <- splinefun(sim_adj_age, f_40yr_percent_d)
plot(1:2); dev.off()
curve(f_40yr_per_d_spline, col='red', from=0, to=82, xlab="years past 40", ylab="percent chance of death")
points(sim_adj_age, f_40yr_percent_d)

# Clamped at infinite rate via pmin
f_40yr_drate <- function(t) inst_rate(pmin(f_40yr_per_d_spline(t), 1),1)
curve(f_40yr_drate, from=0, to=90)

####################################
# Integrations of death rates for exposure calculations in delay
# Now, a special function used in delay equation (Had to put upper bound at 81)
# F_40yr_drate_5yr <- Vectorize(function(t)
# {
#   integrate(f_40yr_drate, lower=max(t-5, 0), upper=min(t, 81))$value
# })
# 
# 
# F_40yr_drate_1yr <- Vectorize(function(t)
# {
#   integrate(f_40yr_drate, lower=max(t-1, 0), upper=min(t, 81))$value
# })


# Does a spline work faster?

# x <- 0:160 / 2
# y <- exp(-5*params['r_b'] - F_40yr_drate_5yr(x))
# plot(x, y, typ='l')
# f <- splinefun(x, y)
# curve(f, add=TRUE, col='red', lty=2)
# F_40yr_drate_5yr <- f
# 
# 
# x <- 0:160 / 2
# y <- exp(-F_40yr_drate_1yr(x))
# plot(x, y, typ='l')
# f <- splinefun(x, y)
# curve(f, add=TRUE, col='red', lty=2)
# F_40yr_drate_1yr <- f

# This is for doing numberical integration of a set of numbers at an even interval
alt_simp_coef <- function(i)
{
  if (i < 8) stop("Invalid Simpson coefficient size")
  
  # pg.117 4.1.14, Numerical Recipes in C, 1st edition
  c(17/48, 50/48, 43/48, 49/48, rep(1, i-8), 49/48, 43/48, 50/48, 17/48) 
}

###################################
# Numerical Delay Differential Equation
Multi <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    
    # # Use table for death_prob, Female 40 (offset 1)
    # 
    # if(is.infinite(r_d)) r_d <- 1e16 # A really large number
    # 
    # # Event B stops at time 5 years after event A (delay equation)
    # dd_b <- if (t < 5 || t> 10) 0 else (1-r_ad)*r_a*lagvalue(t-5, 2)*F_40yr_drate_5yr(t)
    # 
    # # Event A stops at time t=5 years
    # if(t > 5) r_a <- 0
    r_d <- f_40yr_drate(t)
    
    # Interaction is primarily through death (other is panel test)
    liv  <- h_u1+h_t1+a_p1+a_a1+b_n1 # Living (should be same in all models i)
    r_d1 <- r_d + (r_b2*p_bd2*a_p2 +r_b2*rr_b2*p_bd2*a_a2)/liv
    r_d2 <- r_d + (r_b1*p_bd1*a_p1 +r_b1*rr_b1*p_bd1*a_a1)/liv
    r_p1 <- p_p*p_o*r_a2*h_u2 / liv
    r_p2 <- p_p*p_o*r_a1*h_u1 / liv
    
    list(c(
      
      h_u1  = -r_p1*h_u1 -r_a1*h_u1                                                                                    -r_d1*h_u1,
      h_t1  = +r_p1*h_u1 -r_a1*h_t1                                                                                    -r_d1*h_t1,
      a_p1  =            +r_a1*(1-p_o*p_g1)*h_u1+r_a1*(1-p_g1*p_r)*h_t1-r_b1*a_p1                                      -r_d1*a_p1,
      a_a1  =            +r_a1*p_o*p_g1*h_u1    +r_a1*p_g1*p_r*h_t1                         -r_b1*rr_b1*a_a1           -r_d1*a_a1,
      b_n1  =                                                          +r_b1*(1-p_bd1)*a_p1 +r_b1*rr_b1*(1-p_bd1)*a_a1 -r_d1*b_n1,
      b_d1  =                                                          +r_b1*p_bd1*a_p1     +r_b1*rr_b1*p_bd1*a_a1,
      a_c1  =            +r_a1*h_u1             +r_a1*h_t1,
      b_c1  =                                                          +r_b1*a_p1           +r_b1*rr_b1*a_a1,
        
      h_u2  = -r_p2*h_u2 -r_a2*h_u2                                                                                    -r_d2*h_u1,
      h_t2  = +r_p2*h_u2 -r_a2*h_t2                                                                                    -r_d2*h_t1,
      a_p2  =            +r_a2*(1-p_o*p_g2)*h_u2+r_a2*(1-p_g2*p_r)*h_t2-r_b2*a_p2                                      -r_d2*a_p1,
      a_a2  =            +r_a2*p_o*p_g2*h_u2    +r_a2*p_g2*p_r*h_t2                         -r_b2*rr_b2*a_a2           -r_d2*a_a1,
      b_n2  =                                                          +r_b2*(1-p_bd2)*a_p2 +r_b2*rr_b2*(1-p_bd2)*a_a2 -r_d2*b_n1,
      b_d2  =                                                          +r_b2*p_bd2*a_p2     +r_b2*rr_b2*p_bd2*a_a2,
      a_c2  =            +r_a2*h_u2             +r_a2*h_t2,
      b_c2  =                                                          +r_b2*a_p2           +r_b2*rr_b2*a_a2,
      
      disc  = -disc_rate*disc            # Simple discount rate
    ))
  })
}

yinit <- rep(0, 17)
names(yinit) <- c("h_u1", "h_t1", "a_p1", "a_a1", "b_n1", "b_d1", "a_c1", "b_c1", "h_u2", "h_t2", "a_p2", "a_a2", "b_n2", "b_d2", "a_c2", "b_c2", "disc")
yinit[1]  <- 1
yinit[9]  <- 1
yinit[17] <- 1

times <- seq(0, 40, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- ode(yinit, times, Multi, params)))

plot(out)

l1 <- out[,'h_u1'] + out[,'h_t1'] + out[,'a_p1'] + out[,'a_a1'] + out[,'b_n1']
l2 <- out[,'h_u2'] + out[,'h_t2'] + out[,'a_p2'] + out[,'a_a2'] + out[,'b_n2']

plot(l1-l2, typ='l', main="Abs Total Numerical Error")

cat("Consistent population?", all(abs(l1-l2)<1e-8), '\n')

stop("Working Halt")

costs <- function(solution, params)
{
  n <- length(solution[,1])
  simpson <- alt_simp_coef(n)

  with(as.list(params), {
    # Compute Discounted Cost
    cost <- c_a*sum(diff(solution[,'a'])*solution[2:n,'disc']) + # Cost * Number of events in bucket a
            c_b*sum(diff(solution[,'b'])*solution[2:n,'disc']) + # Cost * Number of events in bucket b
            c_t*solution[1, 'h']  # Therapy Cost * Initial healthy individuals
    
    # Step size of simulation
    step     <- solution[2,'time'] - solution[1,'time']
    
    # Total possible life units is integral of discounted time
    life <- sum(simpson*solution[,'disc'])*step
    
    # Permanent disutility for A (integration)
    disA <- d_a*sum(simpson*solution[,'e10']*solution[,'disc'])*step + 
            d_a*sum(simpson*solution[,'e15']*solution[,'disc'])*step +
            d_a*sum(simpson*solution[,'e2' ]*solution[,'disc'])*step
    
    # Event B
    disB <- d_b*sum(simpson*solution[,'db']*solution[,'disc'])*step

    # Death disutility
    disD <- sum(simpson*solution[,'d']*solution[,'disc'])*step       # Death disutility (integration)
    
    c(cost       = unname(cost),
      qaly       = unname(life-disA-disB-disD),
      possible   = unname(life),
      disutility = unname(disA+disB+disD),
      a_count    = unname(solution[n,'a']),
      disutil_a  = unname(disA),
      b_count    = unname(solution[n,'b']),
      disutil_b  = unname(disB),
      dead_count = unname(solution[n,'d']), 
      disutil_d  = unname(disD),
      living     = unname(solution[n,'h']+solution[n,'e10']+solution[n,'e15']+solution[n,'e2'])
      )
  })
}

expected <- function(params) costs(dede(yinit, times, Simple, params), params)

round(expected(params), 4)

params['r_a'] <- params['r_a']*0.5
params['c_t']  <- 2000
round(expected(params), 4)

# Some diagnostic plots
# 
# dev.off()
# plot(out[,'time'],out[,'db'], typ='l', xlim=c(0, 12))
# abline(h=0, col='red', lty=2)
# abline(v=11, col='blue', lty=2)
# abline(v=5, col='green', lty=2)
# 
# dev.off()
# plot(out[,'time'],out[,'b'], typ='l', xlim=c(0, 12))
# abline(v=10, col='green', lty=2)
# 
# dev.off()
# plot(out[,'time'],out[,'e10'], typ='l', xlim=c(0, 12))
# abline(v=10, col='green', lty=2)