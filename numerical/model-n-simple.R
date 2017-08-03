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
      r_p <- 0
      for(j in 1:n)
      {
        if(j != i)
        {
          r_d <- r_d + (r_b[j]*p_bd[j]*(y[map("a_p", j)]+rr_b[j]*y[map("a_a", j)]))/liv
          r_p <- r_p + p_p*p_o*r_a[j]*y[map("h_u", j)] / liv
        }
      }

      rate <- c(rate, 
        (-r_p-r_a[i]-r_d)*y[map("h_u", i)],
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

yinit <- rep(0, length(key)*params[["n"]]+1)
yinit <- c(yinit, 1)
for(i in 1:params[["n"]]) yinit[map("h_u", i)] <- 1


times <- seq(0, 40, by=1/365)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- dede(yinit, times, genModel, params)))

plot(out)

costs <- function(solution, params)
{
  k        <- length(solution[,1])
  simpson  <- alt_simp_coef(k)
  step     <- solution[2,'time'] - solution[1,'time']
  
  with(as.list(params), {
    disc  <- length(key) * n + 3

    cost <- c_t*sum(diff(solution[,disc-1])*solution[2:k,disc]) # Testing costs
    for(i in 1:n)
    {
      # Compute Discounted Cost
      cost <- cost + 
              c_a[i]  *sum(diff(solution[,maps('a_c',i)])*solution[2:k,disc]) +
              c_bs[i] *sum(diff(solution[,maps('b_c',i)])*solution[2:k,disc]) +
              c_bd[i] *sum(diff(solution[,maps('b_d',i)])*solution[2:k,disc]) +
              c_tx[i] *365*sum(simpson*solution[,maps('a_p',i)]*solution[,disc])*step +
              c_alt[i]*365*sum(simpson*solution[,maps('a_a',i)]*solution[,disc])*step +
              c_tx[i] *365*sum(simpson*solution[,maps('b_p',i)]*solution[,disc])*step +
              c_alt[i]*365*sum(simpson*solution[,maps('b_a',i)]*solution[,disc])*step
    }
    
    # Total living in model
    life <- solution[,maps('h_u',i)] +
            solution[,maps('h_t',i)] +
            solution[,maps('a_p',i)] +
            solution[,maps('a_a',i)] +
            solution[,maps('b_p',i)] + 
            solution[,maps('b_a',i)] 

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

    c(cost       = unname(cost),
      qaly       = unname(pQALY-disA-disB),
      possible   = unname(pQALY),
      disutil_a  = unname(disA),
      disutil_b  = unname(disB),
      living     = unname(life[k])
      )
  })
}

expected <- function(params) costs(dede(yinit, times, genModel, params), params)

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