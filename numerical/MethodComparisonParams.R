source("model-n-simple.R")

# Gompertz model for 40yr female
shape <- 0.1007511
rate  <- 0.0008370717

# Using exact gompertz approximation from DES

f_40yr_drate <- function(t) rate*exp(shape*t)

inst_rate <- function(percent, timeframe)
{
  - log(1-percent) / timeframe
}

times <- seq(0, 80, by=1/365)

params       <- list()

params$n    = 1   # Only 1 drug under consideration
params$p_p  = 0.0 # Not reactive panel
params$p_o  = 1.0 # Probability of ordering test
params$p_r  = 1.0 # Probability pre-existing test is read
params$p_bd = 0.05 # Probability of death from B
params$p_g  = 0.2 # Probability of genetic variant
params$r_a  = inst_rate(0.1, 10) # 10% Rate of A over a 10 year period
params$r_b  = inst_rate(0.02, 1) # 2% Rate of B
params$rr_b = 0.7 # Reduced relative risk of B

# Costs
params$c_a   <- 10000 # Cost of Event A
params$c_bs  <- 25000 # Cost of Event B survival
params$c_bd  <- 15000 # Cost of Event B death
params$c_tx  <- 0.5   # Cost of normal treatment
params$c_alt <- 5     # Cost of alternate treatment
  
params$c_t   <- 0   # Cost of test

params$disc_rate <- 0.03 # Discount rate                 

params$d_a   <- 0.05     # Disutility of A
params$d_at  <- 1        # Duration of A in years.
params$d_b   <- 0.1      # Disutility of B
  
params$p_o <- 0.0 # Standard
init    <- generate.initial("none", params)
x       <- dede(init, times, genModel, params)
standard <- costs(x, params)

params$p_o <- 1.0 # Genotype
init    <- generate.initial("single", params)
x       <- dede(init, times, genModel, params)
genotype <- costs(x, params)

print(round(standard,5))
print(round(genotype,5))
print(round( (standard[1] - genotype[1]) / (standard[2] - genotype[2]) , 1))

