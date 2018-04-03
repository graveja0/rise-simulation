library(randtoolbox)

inst_rate <- function(percent, timeframe)
{
  - log(1-percent) / timeframe
}

source("model-n-simple.R")

# Gompertz model for 40yr female
shape <- 0.1007511
rate  <- 0.0008370717

# Using exact gompertz approximation from DES

f_40yr_drate <- function(t) rate*exp(shape*t)

simulate <- function(params)
{

  times <- seq(0, params$horizon, by=1/365)

  params$p_o <- 0.0 # Standard
  init    <- generate.initial("none", params)
  x       <- dede(init, times, genModel, params, control=list(mxhist=1e6))
  standard <- costs(x, params)

  params$p_o <- 1.0 # Genotype
  init    <- generate.initial("single", params)
  x       <- dede(init, times, genModel, params, control=list(mxhist=1e6))
  genotype <- costs(x, params)
  
  x <- c(standard[1], standard[2], genotype[1], genotype[2], (standard[1] - genotype[1]) / (standard[2] - genotype[2]))
  names(x) <- c("std_dCOST", "std_dQALY", "gen_dCOST", "gen_dQALY", "ICER")
  x
}


halton_run <- function(n, start=1)
{
  if(start > 1) halton(start-1, 15, init=TRUE)
  
  apply(halton(n, 15, init=FALSE), 1, function(x) 
  {
    cat(".")
    params       <- list()

    params$n    = 1   # Only 1 drug under consideration
    params$p_p  = 0.0 # Not reactive panel
    params$p_o  = 1.0 # Probability of ordering test
    params$p_r  = 1.0 # Probability pre-existing test is read
    
    
    params$p_bd = x[1] # Probability of death from B
    params$p_g  = x[2] # Probability of genetic variant
    params$r_a  = inst_rate(x[3], 10) # 10% Rate of A over a 10 year period
    params$r_b  = inst_rate(x[4], 1) # 2% Rate of B
    params$rr_b = x[5] # Reduced relative risk of B

    # Costs
    params$c_a   <- exp(4.60517 * x[6] + 6.907755) # Cost of Event A 1000 -> 10000 log scale
    params$c_bs  <- exp(4.60517 * x[7] + 6.907755) # Cost of Event B survival
    params$c_bd  <- exp(4.60517 * x[8] + 6.907755) # Cost of Event B death
    params$c_tx  <- 9.99*x[9]+0.01   # Cost of daily normal treatment
    params$c_alt <- 99.99*x[10]+0.01     # Cost of alternate treatment
    
    params$c_t   <- 1000*x[11]   # Cost of test

    params$disc_rate <- 0.03 # Discount rate                 

    params$d_a   <- 0.5*x[12]     # Disutility of A
    params$d_at  <- 10*x[13]      # Duration of A in years.
    params$d_b   <- 0.5*x[14]      # Disutility of B 
    
    params$horizon <- x[15]*79 + 1 # 1:80 years for simulation
    
    results <- simulate(params)

    x <- with(params, c(p_bd, p_g, r_a, r_b, rr_b, c_a, c_bs, c_bd, c_tx, c_alt, c_t, d_a, d_at, d_b, horizon, results))
    
    names(x) <- c("p_bd", "p_g", "r_a", "r_b", "rr_b", "c_a", "c_bs", "c_bd", "c_tx", "c_alt", "c_t", "d_a", "d_at", "d_b",
                  "horizon", "std_dCOST", "std_dQALY", "gen_dCOST", "gen_dQALY", "ICER")
    x
  })
}

x <- t(halton_run(100))
write.csv(x, "simple-psa-0.csv", row.names=FALSE)

x <- t(halton_run(900, 101))
write.csv(x, "simple-psa-1.csv", row.names=FALSE)

# cat simple-psa-?.csv > psa.csv

x <- read.csv("psa.csv")
x$norm <- (x$ICER - median(x$ICER)) / median(x$ICER)
m <- lm(norm ~ p_bd+p_g+r_a+r_b+rr_b+c_a+c_bs+c_bd+c_tx+c_alt+c_t+d_a+d_at+d_b+horizon, x)

