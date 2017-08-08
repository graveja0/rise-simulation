# Parameters for a 3 conditions model

###################################
# Main model parameters
params <- list(
  n     = 3,
  # Probabilities
  p_o   = 0.75,               # Probability of ordering test
  p_r   = 0.75,               # Probability of reading test if it exists
  p_p   = 1,                  # Probability of panel test when test ordered
  
  p_bd  = c(0.05, 0.05, 0.05),  # Probability of death as direct result of B
  p_g   = c(0.2,  0.2,  0.2),   # Prevalence of targeted gene variant
  
  
  # Risks
  # Rate of healthy having event A over 10 years
  r_a    = c(inst_rate(0.7, 10), inst_rate(0.2, 10), inst_rate(0.1, 10)),
  # Rate of post-A  having event B over 3 years
  r_b    = c(inst_rate(0.2, 3),  inst_rate(0.5, 10), inst_rate(0.5, 10)),
  # Relative Risk for Event B with better treatment informed by genotyping
  rr_b   = c(0.5, 0.5, 0.5),              
  
  # Costs
  c_a    = c(10000, 10000, 10000),# Cost of event A
  c_bs   = c(25000, 25000, 25000),# Cost for event B survival
  c_bd   = c(15000, 15000, 15000),# Cost for event B death
  c_tx   = c(100, 100, 100)/30,   # Daily cost of standard treatment
  c_alt  = c(300, 300, 300)/30,   # Daily cost of alternate treatment
  
  c_t    = 100,                   # Cost of test
  
  # Disutility
  d_a    = c(0.2,  0.2,  0.2),    # Disutility for event A
  d_at   = c(3,    3,    3),      # Time in years for disutility of event A
  d_b    = c(0.15, 0.15, 0.15),   # Permanent Disutility for event B
  
  disc_rate = 1e-12           # For computing discount
)
