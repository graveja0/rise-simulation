source("model-n-simple.R")
load("simulation-results-raw.Rdata")
load("drawn-parameter-values-behavioral-global.Rdata")

y <- subset(results, preemptive == "None" & reactive == "None")

config   <- drawn.parameter.values
times    <- seq(0, 80, by=2/365)
scenario <- "none"
i        <- 1

params   <- generate.params(config, i, scenario)
init     <- generate.initial(scenario, params)
x        <- dede(init, times, genModel, params)

round(costs(x, params), 3)
round(costs(x, params)*2e5, 2)

par(mfrow=c(2,2))
# Incidence of A_SC_A compared
z <- subset(y, resource == "A_SC_A")
n <- length(z$start_time)
plot(z$start_time/365, 1:n/2e5, typ="l", main="Event A1 Compared", ylab="population frac", xlab="years")
lines(x[,1], x[,maps('a_c',1)], col='red', lty=2)
legend(30, 0.2, c("Delay Diff Eq", "Simulation"), lty=c(2,1), col=c("red", "black"))

# Incidence of A_SC_B compared
z <- subset(y, resource == "A_SC_B")
n <- length(z$start_time)
plot(z$start_time/365, (1:n)/2e5, typ="l", main="Event A2 Compared", ylab="population frac", xlab="years")
lines(x[,1], x[,maps('a_c',2)], col='red', lty=2)


# Incidence of B_SC_A compared
z <- subset(y, resource == "B_SC_A")
n <- length(z$start_time)
plot(z$start_time/365, 1:n/2e5, typ="l", main="Event B1 Compared", ylab="population frac", xlab="years")
lines(x[,1], (x[,maps('b_c',1)]+x[,maps('b_d',1)]) , col='red', lty=2)


# Incidence of B_SC_B compared
z <- subset(y, resource == "B_SC_B")
n <- length(z$start_time)
plot(z$start_time/365, 1:n/2e5, typ="l", main="Event B2 Compared", ylab="population frac", xlab="years")
lines(x[,1], (x[,maps('b_c',2)]+x[,maps('b_d',2)] ) , col='red', lty=2)


dev.off()
# B Death Rates compared
z <- subset(y, resource %in% c("B_Death_SC_A", "B_Death_SC_B"))
plot(ecdf(z$start_time/365), main="Event B Death Compared", xlab="cumulative", ylab="years")
lines(x[,1], (x[,maps('b_d',1)]+x[,maps('b_d',2)] ) * 2e5 / length(unique(z$name)) , col='red', lty=2)
