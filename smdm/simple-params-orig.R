inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

##############################################################
##
## Standardized Parameter List for Simple Model Comparison
##
params <- list(
  # Controls for model execution
  n          = 1000000,      # DES simulations to perform
  resolution = 7/365,        # Diff Eq Time step for DEQ approach
  interval   = 1,            # Markov Interval
  horizon    = 40,           # Time horizon of simulation
  wtp        = 100000,       # Willingness to pay threshold

  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = 0.1007511,
  rate    = 0.0008370717,
  
  # Probabilities and rates
  p_o  = 1.0,                # Probability of ordering test
  p_bd = 0.10,               # Probability of death from B
  p_g  = 0.2,                # Probability of genetic variant
  r_a  = inst_rate(0.1, 10), # 10% Rate of A over a 10 year period
  r_b  = inst_rate(0.02, 1), # 2% Rate of B over a 1 year period
  rr_b = 0.7,                # Reduced relative risk of B
  
  # Costs
  c_a   = 10000,             # Cost of Event A
  c_bs  = 25000,             # Cost of Event B survival
  c_bd  = 15000,             # Cost of Event B death
  c_tx  = 0.5,               # Cost of normal treatment
  c_alt = 5,                 # Cost of alternate treatment
  c_t   = 100,               # Cost of test
  
  # Disutilities
  d_a   = 0.05,              # Disutility of A
  d_at  = 1,                 # Duration of A in years.
  d_b   = 0.1,               # Disutility of B
  
  # Discounting
  disc  = 0.03               # Annual Discount Rate
)