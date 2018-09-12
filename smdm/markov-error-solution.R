params <- c(
  x_rate = 0.1,
  y_rate = 0.5,
  lag    = 1
)

out <- dede(c(i=1, x=0, y=0, z=0), times, Truth, params)
x <- out[366:14601, 'x']
y <- times[366:14601] - 1
error <- function(r) sum((x[1]*exp(-r*y) - x)^2)
optimize(function(x) error(x), interval=c(1e-6, 10))$minimum

plot(y, x, typ='l')
lines(y, x[1]*exp(-params$x_rate*y), col='red')

params$x_rate <- 0.2
out <- dede(c(i=1, x=0, y=0, z=0), times, Truth, params)
x <- out[366:14601, 'x']
y <- times[366:14601] - 1
error <- function(r) sum((x[1]*exp(-r*y) - x)^2)
optimize(function(x) error(x), interval=c(1e-6, 10))$minimum
plot(y, x, typ='l')
lines(y, x[1]*exp(-params$x_rate*y), col='red')

params$x_rate <- 0.5
out <- dede(c(i=1, x=0, y=0, z=0), times, Truth, params)
x <- out[366:14601, 'x']
y <- times[366:14601] - 1
error <- function(r) sum((x[1]*exp(-r*y) - x)^2)
optimize(function(x) error(x), interval=c(1e-6, 10))$minimum
plot(y, x, typ='l')
lines(y, x[1]*exp(-params$x_rate*y), col='red')

#

a <- params['x_rate']
b <- params['y_rate']
tau <- params['lag']


f1 <- function(t) if(a==b) a*exp(-a*t)*t else a/(b-a)*(exp(-a*t) - exp(-b*t))

f  <- function(t) a*exp(-a*t) - exp(-b*tau)*a*exp(-a*(t-tau))
s  <- Vectorize(function(t) integrate(function(x) exp(b*x)*f(x), lower=0, upper=t)$value)

f2 <- function(t, c) exp(-b*t)*(c+s(t))

c <- (f1(tau) - f2(tau, 0)) / exp(-b*tau)

## This works
f2 <- function(t, c) exp(-b*t) * (c+
      a/(b-a)*exp((b-a)*t) -
      a/(b-a)*exp((b-a)*(t-tau))
    )
c <- (f1(tau) - f2(tau, 0)) / exp(-b*tau)

## Simplification
f2 <- function(t, c) exp(-b*t) * (c+
      a/(b-a)*exp((b-a)*t)*(1 - exp((a-b)*tau))
    )
c <- (f1(tau) - f2(tau, 0)) / exp(-b*tau)

## Deeper Simplification
f2 <- function(t, c) 
  c*exp(-b*t) + 
  a/(b-a)*(1-exp( (a-b) *tau)) * exp(-a*t)
c <- (f1(tau) - f2(tau, 0)) / exp(-b*tau)
  

X = function(t) ifelse(t <= tau, f1(t), f2(t, c))

curve(X, from=0, to=40, n=365*40)
lines(out[,'time'], out[,'x'], col='red', lty=2)

# WOW This show the core issue. 
# Loading still has risk present.
markov <- function(t) ifelse(t < tau, 0, a*((1-a)^floor(t-tau)))
curve(markov, from=0, to=40)
curve(X, n=365*40, add=TRUE, col='red', lty=2)

corrected <- function(t) ifelse(t < tau, 0, f1(tau)*((1-a)^floor(t-tau)))
curve(corrected, n=365*40, add=TRUE, col='brown', lty=3)


