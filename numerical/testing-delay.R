# This little bit of code is a minimal example of an event with a finite duration
# It creates a delay-differential

library(deSolve)

###################################
# Main model parameters
params <- c(
  alpha = 0.3,
  beta  = 0.01
)

# Incoming Event Rate from healthy population
f <- function(t, params) with(as.list(params), alpha)
# Death function
g <- function(t, params) with(as.list(params), beta*t)
G <- function(t, params) integrate(function(x) g(x, params), lower=max(t-2,0), upper=t)$value


###################################
# Numerical Delay Differential Equation
DelayTest <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    finished <- if(t < 2) 0 else lagderiv(t-2, 2)*exp(-G(t, params))
    
    list(
      c(
        healthy  = -(f(t, params)+g(t, params))*healthy,
        eventA   = f(t, params)*healthy,
        experA   = f(t, params)*healthy -g(t, params)*experA - finished,
        postA    = finished - g(t, params)*postA,
        deaths   = g(t, params)*(healthy+experA+postA)
      )
    )
  })
}

yinit <- c(healthy=1, eventA=0, experA=0, postA=0, deaths=0)

times <- seq(0, 10, by=1/20)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- dede(yinit, times, DelayTest, params)))
plot(100*out)
plot(1 - out[,'healthy'] - out[,'deaths'] - out[,'experA'] - out[,'postA'], typ='l', ylab="Error")