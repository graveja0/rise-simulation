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
  c_a1   = 10000,              # Cost of event A
  c_bs1  = 25000,              # Cost for event B survival
  c_bd1  = 15000,              # Cost for event B death
  c_tx1  = 100/30,             # Daily cost of standard treatment
  c_alt1 = 300/30,             # Daily cots of alternate treatment

  
  c_a2   = 10000,              # Cost of event A
  c_bs2  = 25000,              # Cost for event B survival
  c_bd2  = 15000,              # Cost for event B death
  c_tx2  = 100/30,             # Daily cost of standard treatment
  c_alt2 = 700/30,             # Daily cots of alternate treatment

  c_t    = 100,                # Cost of test
  
  # Disutility
  d_a1   = 0.2,                # Disutility for event A
  d_at1  = 3*365,              # Time in days for disutility of event A
  d_b1   = 0.15,               # Permanent Disutility for event B
  
  d_a2   = 0.2,                # Disutility for event A
  d_at2  = 3*365,              # Time in days for disutility of event A
  d_b2   = 0.15,               # Permanent Disutility for event B
  

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
      b_c1  =                                                          +r_b1*(1-p_bd1)*a_p1 +r_b1*rr_b1*(1-p_bd1)*a_a1,
        
      h_u2  = -r_p2*h_u2 -r_a2*h_u2                                                                                    -r_d2*h_u1,
      h_t2  = +r_p2*h_u2 -r_a2*h_t2                                                                                    -r_d2*h_t1,
      a_p2  =            +r_a2*(1-p_o*p_g2)*h_u2+r_a2*(1-p_g2*p_r)*h_t2-r_b2*a_p2                                      -r_d2*a_p1,
      a_a2  =            +r_a2*p_o*p_g2*h_u2    +r_a2*p_g2*p_r*h_t2                         -r_b2*rr_b2*a_a2           -r_d2*a_a1,
      b_n2  =                                                          +r_b2*(1-p_bd2)*a_p2 +r_b2*rr_b2*(1-p_bd2)*a_a2 -r_d2*b_n1,
      b_d2  =                                                          +r_b2*p_bd2*a_p2     +r_b2*rr_b2*p_bd2*a_a2,
      a_c2  =            +r_a2*h_u2             +r_a2*h_t2,
      b_c2  =                                                          +r_b2*a_p2           +r_b2*rr_b2*a_a2,
      
      test  = +p_o*(r_a1*h_u1 + r_a2*h_u2),
      disc  = -disc_rate*disc            # Simple discount rate
    ))
  })
}

yinit <- rep(0, 18)
names(yinit) <- c("h_u1", "h_t1", "a_p1", "a_a1", "b_n1", "b_d1", "a_c1", "b_c1",
                  "h_u2", "h_t2", "a_p2", "a_a2", "b_n2", "b_d2", "a_c2", "b_c2",
                  "test", "disc")
yinit[1]  <- 1
yinit[9]  <- 1
yinit["disc"] <- 1

times <- seq(0, 40, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- ode(yinit, times, Multi, params)))

plot(out)

# Check that living population totals are in agreement within numerical error
l1 <- out[,'h_u1'] + out[,'h_t1'] + out[,'a_p1'] + out[,'a_a1'] + out[,'b_n1']
l2 <- out[,'h_u2'] + out[,'h_t2'] + out[,'a_p2'] + out[,'a_a2'] + out[,'b_n2']
#plot(l1-l2, typ='l', main="Abs Total Numerical Error")
cat("Consistent population?", all(abs(l1-l2)<1e-8), '\n')



costs <- function(solution, params)
{
  n <- length(solution[,1])
  simpson <- alt_simp_coef(n)
  step     <- solution[2,'time'] - solution[1,'time']
  
  with(as.list(params), {
    # Compute Discounted Cost
    cost <- c_a1  *sum(diff(solution[,'a_c1'])*solution[2:n,'disc']) +
            c_bs1 *sum(diff(solution[,'b_c1'])*solution[2:n,'disc']) +
            c_bd1 *sum(diff(solution[,'b_d1'])*solution[2:n,'disc']) +
            c_tx1 *365*sum(simpson*solution[,'a_p1']*solution[,'disc'])*step +
            c_alt1*365*sum(simpson*solution[,'a_a1']*solution[,'disc'])*step +
            c_a2  *sum(diff(solution[,'a_c2'])*solution[2:n,'disc']) +
            c_bs2 *sum(diff(solution[,'b_c2'])*solution[2:n,'disc']) +
            c_bd2 *sum(diff(solution[,'b_d2'])*solution[2:n,'disc']) +
            c_tx2 *365*sum(simpson*solution[,'a_p2']*solution[,'disc'])*step +
            c_alt2*365*sum(simpson*solution[,'a_a1']*solution[,'disc'])*step +
            c_t   *sum(diff(solution[,'test'])*solution[2:n,'disc'])

    # Total possible life units is integral of discounted time
    life <- sum(simpson*solution[,'disc'])*step
    
    disA <- NA # Need to figure out
    
    # Permanent disutility for B (integration)
    disB <- d_b1*sum(simpson*solution[,'b_n1']*solution[,'disc'])*step + 
            d_b2*sum(simpson*solution[,'b_n2']*solution[,'disc'])*step
      
    # Total living in model
    life <- solution[,'h_u1'] +
            solution[,'h_t1'] +
            solution[,'a_p1'] +
            solution[,'a_a1'] +
            solution[,'b_n1']
    pQALY <- sum(simpson*life*solution[,'disc'])*step

    c(cost       = unname(cost),
      qaly       = unname(pQALY-disA-disB),
      possible   = unname(pQALY),
      disutil_a  = unname(disA),
      disutil_b  = unname(disB),
      a1_count   = unname(solution[n,'a_c1']),
      b1_count   = unname(solution[n,'b_n1']+solution[n,'b_d1']),
      living     = unname(life[n])
      )
  })
}

expected <- function(params) costs(ode(yinit, times, Multi, params), params)

round(expected(params), 4)
stop("Working Halt")
# params['r_a'] <- params['r_a']*0.5
# params['c_t']  <- 2000
# round(expected(params), 4)

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