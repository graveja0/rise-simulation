library(deSolve)

ss_death <- read.csv("ss-death-2011.csv")

source('params-3.R')

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

###################################
# Numerical Delay Differential Equation
genModel <- function(params)
{
  function(t, y, ignore)
  {
    with(as.list(c(y, params)), {
    
      r_d <- f_40yr_drate(t)
    
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
      
      #### FIXME!!!! 
      r_d <- r_d + (r_b2*p_bd2*(a_p2+rr_b2*a_a2))/liv
      r_p <- p_p*p_o* r_a2*h_u2 / liv
      ####
      
      rate <- c(rate, 
        (-r_p-r_a[i]-r_d1[i])*y[map("h_u", i)],
        r_p*y[map("h_u", i)] + (-r_a[i]-r_d)*y[map("h_t", i)],
        r_a[i]*(1-p_o*p_g[i])*y[map("h_u", i)]+r_a[i]*(1-p_g[i]*p_r)*y[map("h_t", i)]-r_b[i]*y[map("a_p", i)] -r_d*y[map("a_p", i)],
        r_a[i]*p_o*p_g[i]*y[map("h_u", i)]+r_a[i]*p_g[i]*p_r*y[map("h_t", i)] -r_b[i]*rr_b[i]*y[map("a_a", i)] -r_d*y[map("a_a", i)],
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
    
    tests <- p_o*sum(sapply(1:n, function(i) {
      r_a[i]*y[map("h_u", i)]
    }))
      
    list(c(rate,
           tests,
           -disc_rate*y[n * length(key) + 2]
           ))
    })
  }
}

yinit <- rep(0, length(key)*params[["n"]]+1)
yinit <- c(yinit, 1)
for(i in 1:params[["n"]]) yinit[map("h_u1", i)] <- 1


times <- seq(0, 40, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- dede(yinit, times, genModel(params), NULL)))

plot(out)

# Check that living population totals are in agreement within numerical error
l1 <- out[,'h_u1'] + out[,'h_t1'] + out[,'a_p1'] + out[,'a_a1'] + out[,'b_p1'] + out[,'b_a1']
l2 <- out[,'h_u2'] + out[,'h_t2'] + out[,'a_p2'] + out[,'a_a2'] + out[,'b_p2'] + out[,'b_a2']
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
            c_tx1 *365*sum(simpson*solution[,'b_p1']*solution[,'disc'])*step +
            c_alt1*365*sum(simpson*solution[,'b_a1']*solution[,'disc'])*step +
      
            c_a2  *sum(diff(solution[,'a_c2'])*solution[2:n,'disc']) +
            c_bs2 *sum(diff(solution[,'b_c2'])*solution[2:n,'disc']) +
            c_bd2 *sum(diff(solution[,'b_d2'])*solution[2:n,'disc']) +
            c_tx2 *365*sum(simpson*solution[,'a_p2']*solution[,'disc'])*step +
            c_alt2*365*sum(simpson*solution[,'a_a2']*solution[,'disc'])*step +
            c_tx2 *365*sum(simpson*solution[,'b_p2']*solution[,'disc'])*step +
            c_alt2*365*sum(simpson*solution[,'b_a2']*solution[,'disc'])*step +
      
            c_t   *sum(diff(solution[,'test'])*solution[2:n,'disc'])

    
    # Total living in model
    life <- solution[,'h_u1'] +
            solution[,'h_t1'] +
            solution[,'a_p1'] +
            solution[,'a_a1'] +
            solution[,'b_p1'] + 
            solution[,'b_a1'] 

    # Total possible life units is integral of discounted time
    pQALY <- sum(simpson*life*solution[,'disc'])*step

    # Temp disutility of A
    disA <- d_a1*sum(simpson*solution[,'a_q1']*solution[,'disc'])*step + 
            d_a2*sum(simpson*solution[,'a_q2']*solution[,'disc'])*step
      
    # Permanent disutility for B (integration)
    disB <- d_b1*sum(simpson*solution[,'b_p1']*solution[,'disc'])*step + 
            d_b1*sum(simpson*solution[,'b_a1']*solution[,'disc'])*step +
            d_b2*sum(simpson*solution[,'b_p2']*solution[,'disc'])*step +
            d_b2*sum(simpson*solution[,'b_a2']*solution[,'disc'])*step 

    c(cost       = unname(cost),
      qaly       = unname(pQALY-disA-disB),
      possible   = unname(pQALY),
      disutil_a  = unname(disA),
      disutil_b  = unname(disB),
      a1_count   = unname(solution[n,'a_c1']),
      b1_count   = unname(solution[n,'b_c1']+solution[n,'b_d1']),
      living     = unname(life[n])
      )
  })
}

expected <- function(params) costs(dede(yinit, times, Multi, params), params)

sol <- round(expected(params), 4)
print(sol)

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