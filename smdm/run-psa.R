source("simple-deq.R")
source("simple-des.R")
source("simple-markov.R")
source("simple-params.R")

library(randtoolbox)
library(microbenchmark)
library(progress)

s1   <- function(x, total=50)  total*x
s2   <- function(x, total=50)  total - total*x
draw <- function(x, n=1) rbeta(n, s1(x), s2(x))

psa_run <- function(n, FUN=deq_icer, seed=314152)
{
  set.seed(seed)

  pb <- progress_bar$new(format= "(:spin) [:bar] :percent\r", total = n)
  
  result <- sapply(1:n, function(x) 
  {
    pb$tick()
    
    params$p_bd = draw(0.1)  # Probability of death from B
    params$p_g  = draw(0.75) # Probability of genetic variant
    params$r_a  = draw(0.6)  # Inst Rate of A over a 10 year period
    params$r_b  = draw(0.7)  # Inst Rate of B
    params$rr_b = draw(0.1)  # Reduced relative risk of B
    
    # Costs
    params$c_a   <- exp(draw(0.5)-0.5 + log(5800)) # Cost of Event A 1000 -> 10000 log scale
    params$c_bs  <- exp(draw(0.5)-0.5 + log(1200)) # Cost of Event B survival
    params$c_bd  <- exp(draw(0.5)-0.5 + log(27500)) # Cost of Event B death
    params$c_tx  <- exp(draw(0.5)-0.5 + log(4))   # Cost of daily normal treatment
    params$c_alt <- exp(draw(0.5)-0.5 + log(33))     # Cost of alternate treatment
    params$c_t   <- exp(draw(0.5)-0.5 + log(335))   # Cost of test

    params$d_a   <- draw(0.16)    # Disutility of A
    params$d_at  <- 9.5*draw(0.5) + 0.5  # Duration of Disutility A in years.
    params$d_b   <- draw(0.12)     # Disutility of B 
    
    params$horizon <- draw(0.5)*79 + 1 # 1:80 years for simulation
    
    et <- microbenchmark(results <- FUN(params), times=1L)
    
    x <- with(params, c(p_bd, p_g, r_a, r_b, rr_b, c_a, c_bs, c_bd, c_tx, c_alt, c_t, d_a, d_at, d_b, horizon, et$time/1e6, results))
    
    names(x) <- c("p_bd", "p_g", "r_a", "r_b", "rr_b", "c_a", "c_bs", "c_bd", "c_tx", "c_alt", "c_t", "d_a", "d_at", "d_b",
                  "horizon", "system.time", names(results))
    x
  })
  pb$terminate()
  
  result
}

#x <- t(psa_run(5000, deq_icer))
#write.csv(x, "data/deq-psa-middle.csv", row.names=FALSE)
# 
#x <- t(psa_run(5000, des_icer))
#write.csv(x, "data/des-psa-middle.csv", row.names=FALSE)
# 
#x <- t(psa_run(5000, markov_icer))
#write.csv(x, "data/markov-psa-middle.csv", row.names=FALSE)

