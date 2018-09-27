inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

##############################################################
##
## Standardized Parameter List for Simple Model Comparison
##
params <- list(
  # Controls for model execution
  n          = 100000,      # DES simulations to perform
  resolution = 7/365,        # Diff Eq Time step for DEQ approach
  interval   = 1,            # Markov Interval
  horizon    = 40,           # Time horizon of simulation
  wtp        = 100000,       # Willingness to pay threshold

  # Gompertz model of secular death for 40yr female
  # fit from 2012 social security data
  shape   = 0.1007511,
  rate    = 0.0008370717,
  
  # Probabilities and rates
  p_o  = 1.0,                # Probability of ordering test (overwritten by runs to 0 and 1)
  p_bd = 0.1,                # Probability of death from B
  p_g  = 0.75,               # Probability of genetic variant (majority)
  r_a  = 0.6,                # Inst Rate of A
  r_b  = 0.7,                # Inst Rate of B 
  rr_b = 0.1,                # Reduced relative risk of B
  
  # Costs
  c_a   = 5800,              # Cost of Event A
  c_bs  = 1200,              # Cost of Event B survival
  c_bd  = 27500,             # Cost of Event B death
  c_tx  = 4,                 # Cost of normal treatment
  c_alt = 33,                # Cost of alternate treatment
  c_t   = 335,               # Cost of test
  
  # Disutilities
  d_a   = 0.16,               # Disutility of A
  d_at  = 7.5,               # Duration of A in years.
  d_b   = 0.12,              # Disutility of B
  
  # Discounting
  disc  = 0.03               # Annual Discount Rate
)