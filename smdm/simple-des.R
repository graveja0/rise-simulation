library(flexsurv)

inst_rate <- function(percent, timeframe) -log(1-percent) / timeframe

  ##############################################################
 ##
## Standardized Parameter List for Simple Model Comparison
##
params <- list(
  # Control
  n       = 100000,              # DES simulations to perform
  horizon = 40,              # Time horizon of simulation
  
  # Gompertz model of secular death for 40yr female
  shape   = 0.1007511,
  rate    = 0.0008370717,
  
  # Probabilities and rates
  p_o  = 1.0,                # Probability of ordering test
  p_bd = 0.05,               # Probability of death from B
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
  disc  = 0.03           # Annual
)

# Perform all random draws
DES <- with(params, {
  x <- data.frame(
    name          = paste0("patient",1:n),
    variant       = sample(c(TRUE, FALSE), size=n, replace=TRUE, prob=c(p_g, 1-p_g)),
    secular_death = rgompertz(n, shape, rate)*365,
    indication       = rexp(n, r_a / 365)
  )

  x$indication[x$secular_death > horizon*365]    <- NA
  x$secular_death[x$secular_death > horizon*365] <- NA
  x$indication[x$secular_death < x$indication]   <- NA
  x$tested <- !is.na(x$indication) * rbinom(n, 1, p_o)
  x$treat <- NA
  x$treat[!is.na(x$indication)] <- "Primary"
  x$treat[!is.na(x$indication) & x$tested & x$variant] <- "Alternate"
  
  rr <- ifelse(x$treat=="Alternate", rr_b, 1.0)
  x$adverse <- x$indication + suppressWarnings(rexp(n, r_b*rr/365))
  x$adverse[x$adverse > horizon*365] <- NA
  x$adverse[x$secular_death < x$adverse] <- NA

  x$adverse_death <- rbinom(n, 1, p_bd) * x$adverse
  x$adverse_death[is.na(x$adverse) | x$adverse_death == 0] <- NA
  
  x$secular_death[x$secular_death > x$adverse_death] <- NA
  
  x$death <- x$secular_death
  x$death[!is.na(x$adverse_death)] <- x$adverse_death[!is.na(x$adverse_death)]
  
  x$end_of_sim <- x$death
  x$end_of_sim[is.na(x$end_of_sim)] <- horizon*365
  
  x
})

# Discounted integral from A to B
discount_int <- function(ar, A, B)
{
  r <- (1 + ar)^(1/365)-1
  (exp(-r*A)-exp(-r*B)) / r
}

# Discounted integrals summed
disum <- function(disc, A, B) sum(discount_int(disc, A, B), na.rm=TRUE)

# Discounted discrete sum
# Each value is a time of event in units of days (365 per year)
ddsum <- function(x, drate) sum(1 / ((1 + drate)^(x/365)), na.rm=TRUE)

des_summary <- function(df, params)
{
  with(c(df, params), {
    # Testing costs
    test.cost  <- c_t  * ddsum(indication, disc)
    treat.cost <- c_a  * ddsum(indication, disc) + 
                  c_bs * ddsum(adverse[is.na(adverse_death)], disc) +
                  c_bd * ddsum(adverse_death, disc) 

    pri <- !is.na(treat) & treat=="Primary"
    alt <- !is.na(treat) & treat=="Alternate"
    drug.cost <- c_tx *disum(disc, indication[pri], end_of_sim[pri]) +
                 c_alt*disum(disc, indication[alt], end_of_sim[alt])

    # Total living in model
    life <- (n - sum(!is.na(death)))

    # Total possible discounted life units is integral of discounted time
    pQALY <- disum(disc, 0, end_of_sim)
 
    # Temp disutility of Indication
    cutoff <- indication + 365*d_at
    cutoff <- ifelse(!is.na(cutoff) & !is.na(adverse) & adverse < cutoff, adverse, cutoff)
    cutoff <- ifelse(!is.na(cutoff) & !is.na(adverse) & secular_death < cutoff, secular_death, cutoff)
    disA <- d_a*disum(disc, indication, cutoff)
    
    # Permanent disutility for Event
    disB <- d_b*disum(disc, adverse, end_of_sim)

    c(dCOST       = unname(treat.cost+test.cost+drug.cost),
      dQALY       = unname(pQALY-disA-disB),
      possible    = unname(pQALY),
      fatal_b     = sum(!is.na(adverse_death)),
      living      = unname(life),
      disutil_a   = unname(disA),
      disutil_b   = unname(disB),
      dCOST.test  = unname(test.cost),
      dCOST.drug  = unname(drug.cost),
      dCOST.treat = unname(treat.cost)
      )/n
  })
}

head(DES)
des_summary(DES, params)