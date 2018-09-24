source("simple-deq.R")
source("simple-des.R")
source("simple-markov.R")
source("simple-params.R")

library(randtoolbox)
library(microbenchmark)
library(progress)


halton_run <- function(n, FUN=deq_icer, start=1)
{
  if(start > 1) halton(start-1, 15, init=TRUE)
  
  pb <- progress_bar$new(format= "(:spin) [:bar] :percent\r", total = n)
  
  result <- apply(halton(n, 15, init=(start == 1)), 1, function(x) 
  {
    pb$tick()
    
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
               
    params$d_a   <- 0.5*x[12]     # Disutility of A
    params$d_at  <- 10*x[13]      # Duration of A in years.
    params$d_b   <- 0.5*x[14]     # Disutility of B 
    
    params$horizon <- x[15]*79 + 1 # 1:80 years for simulation
    
    #cat("running ")
    #print(params)
    #cat("\n")
    et <- microbenchmark(results <- FUN(params), times=1L)
    
    x <- with(params, c(p_bd, p_g, r_a, r_b, rr_b, c_a, c_bs, c_bd, c_tx, c_alt, c_t, d_a, d_at, d_b, horizon, et$time/1e6, results))
    
    names(x) <- c("p_bd", "p_g", "r_a", "r_b", "rr_b", "c_a", "c_bs", "c_bd", "c_tx", "c_alt", "c_t", "d_a", "d_at", "d_b",
                  "horizon", "system.time", names(results))
    x
  })
 pb$terminate()
  
  result
}

#x <- t(halton_run(5000, deq_icer))
#write.csv(x, "deq-icer-psa-0.csv", row.names=FALSE)

#x <- t(halton_run(5000, des_icer))
#write.csv(x, "des-icer-psa-0.csv", row.names=FALSE)

x <- t(halton_run(5000, markov_icer))
write.csv(x, "des-markov-psa-0.csv", row.names=FALSE)


# write.csv(x, "simple-psa-0.csv", row.names=FALSE)
# 
# x <- t(halton_run(900, 101))
# write.csv(x, "simple-psa-1.csv", row.names=FALSE)

# cat simple-psa-?.csv > psa.csv

# x <- read.csv("psa.csv")
# x$norm <- (x$ICER - median(x$ICER)) / median(x$ICER)
# m <- lm(norm ~ p_bd+p_g+r_a+r_b+rr_b+c_a+c_bs+c_bd+c_tx+c_alt+c_t+d_a+d_at+d_b+horizon, x)

