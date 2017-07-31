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

###################################
# Numerical Delay Differential Equation
DelayTest <- function(t, y, params)
{
  with(as.list(c(y, params)), {
    finished <- if(t < 2) 0 else lagderiv(t-2, 2)*exp(-exposure)
    
    list(
      c(
        healthy  = -(f(t, params)+g(t, params))*healthy,
        eventA   = f(t, params)*healthy,
        experA   = f(t, params)*healthy -g(t, params)*experA - finished,
        postA    = finished - g(t, params)*postA,
        deaths   = g(t, params)*(healthy+experA+postA),
        
        # Alternate method of determining integral G when unavailable
        accloss  = g(t, params),  # Accumulator
        exposure = g(t, params) - if(t < 2) 0 else lagderiv(t-2, 6) # Lagged Interval
      )
    )
  })
}

yinit <- c(healthy=1, eventA=0, experA=0, postA=0, deaths=0, accloss=0, exposure=0)

times <- seq(0, 10, by=1/20)  # units of years, increments of days, everyone dies after 120, so simulation is cut short
print(system.time(out <- dede(yinit, times, DelayTest, params)))

plot(100*out)
plot(1 - out[,'healthy'] - out[,'deaths'] - out[,'experA'] - out[,'postA'], typ='l', ylab="Error")

pp <- function(data, vars)
{
  colors <- rainbow_hcl(length(vars))
  plot(data[,'time'], data[,vars[1]], typ='l', ylab='', xlab="time", col=colors[1])
  for(i in 2:length(vars))
  {
    lines(data[,'time'], data[,vars[i]], col=colors[i])
  }
}

#pp(100*out, c("altG", "G"))