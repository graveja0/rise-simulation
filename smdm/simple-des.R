library(flexsurv)

#source("simple-params.R")

# Perform all random draws in the manner of Discrete Event Simulation
# This can be done without larger frameworks due to simplicity of model
des_simulation <- function(params)
{
  with(params, {
    x <- data.frame(
      name          = paste0("patient",1:n),
      variant       = sample(c(TRUE, FALSE), size=n, replace=TRUE, prob=c(p_g, 1-p_g)),
      secular_death = rgompertz(n, shape, rate)*365,
      indication    = rexp(n, r_a / 365)
    )
    

    x$secular_death[x$secular_death > horizon*365] <- NA
    x$indication[x$indication > horizon*365]       <- NA
    x$indication[x$secular_death < x$indication]   <- NA
    x$tested <- !is.na(x$indication)  & rbinom(n, 1, p_o)
    x$treat  <- NA
    x$treat[!is.na(x$indication)] <- "Primary"
    x$treat[!is.na(x$indication) & x$tested & x$variant] <- "Alternate"
    
    rr <- ifelse(x$treat=="Alternate", rr_b, 1.0)
    x$adverse <- x$indication + suppressWarnings(rexp(n, r_b*rr/365))
    x$adverse[is.nan(x$adverse)] <- NA
    x$adverse[x$adverse > horizon*365] <- NA
    x$adverse[x$secular_death < x$adverse] <- NA
    
    x$adverse_death <- rbinom(n, 1, p_bd) * x$adverse
    x$adverse_death[is.na(x$adverse) | x$adverse_death == 0] <- NA
    
    x$secular_death[x$secular_death > x$adverse_death] <- NA
    
    x$death <- x$secular_death
    x$death[!is.na(x$adverse_death)] <- x$adverse_death[!is.na(x$adverse_death)]
    
    x$end_of_sim <- x$death
    x$end_of_sim[is.na(x$end_of_sim)] <- horizon*365
    
    # Useful for computing temp disutility of A
    x$cutoff <- x$indication + 365*d_at
    x$cutoff <- ifelse(!is.na(x$cutoff) & !is.na(x$adverse) & x$adverse < x$cutoff, x$adverse, x$cutoff)
    x$cutoff <- ifelse(!is.na(x$cutoff) & !is.na(x$secular_death) & x$secular_death < x$cutoff, x$secular_death, x$cutoff)
    x$cutoff[x$cutoff > horizon*365] <- horizon*365
    
    x
  })
}

# Discounted integral from A to B at annual rate ar of discounting
discount_int <- function(ar, A, B)
{
  r <- (1 + ar)^(1/365)-1
  (exp(-r*A)-exp(-r*B)) / r
}

# Discounted Integral over interval
disum <- function(disc, A, B) sum(discount_int(disc, A, B), na.rm=TRUE)

# Discounted Discrete sum (aka single event discounting)
# Each value is a time of event in units of days (365 per year)
ddsum <- function(x, drate) sum(1 / ((1 + drate)^(x/365)), na.rm=TRUE)

des_summary <- function(df, params)
{
  with(c(df, params), {
    # Testing costs
    test.cost  <- c_t  * ddsum(indication[tested], disc)
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
    pQALY <- disum(disc, 0, end_of_sim) / 365
 
    # Temp disutility of Indication
    disA   <- d_a*disum(disc, indication, cutoff) / 365
    
    # Permanent disutility for Event
    disB   <- d_b*disum(disc, adverse, end_of_sim) / 365

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

des_icer <- function(params)
{
  params$p_o <- 0.0 # No testing, reference
  reference  <- des_summary(des_simulation(params), params)

  params$p_o <- 1.0 # Genotype testing upon indication
  genotype   <- des_summary(des_simulation(params), params)

  c( ICER       = unname((genotype['dCOST'] - reference['dCOST']) / (genotype['dQALY'] - reference['dQALY'])),
     NMB        = unname((reference['dCOST'] - genotype['dCOST']) + params$wtp*(genotype['dQALY'] - reference['dQALY'])),
     dCOST.ref  = unname(reference['dCOST']),
     dCOST.test = unname(genotype['dCOST']),
     dQALY.ref  = unname(reference['dQALY']),
     dQALY.test = unname(genotype['dQALY'])
  )
}
