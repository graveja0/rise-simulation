setwd("numerical")

load("drawn-parameter-values-behavioral-global.Rdata")

source("model-n-simple.R")

config   <- drawn.parameter.values
times    <- seq(0, 40, by=1/365)
scenario <- "none"

m <- model.run(config, 1, "none")
d <- data.frame(dQALY=m['dQALY'], dCOST=m['dCOST'], stategy="None", fatal_b = 2e+5*m['fatal_b'])

m <- model.run(config, 1, "reactive-single")
d <- rbind(d, data.frame(dQALY=m['dQALY'], dCOST=m['dCOST'], stategy="Reactive Single", fatal_b = 2e+5*m['fatal_b']))

m <- model.run(config, 1, "reactive-panel")
d <- rbind(d, data.frame(dQALY=m['dQALY'], dCOST=m['dCOST'], stategy="Reactive Panel", fatal_b = 2e+5*m['fatal_b']))

m <- model.run(config, 1, "preemptive-panel")
d <- rbind(d, data.frame(dQALY=m['dQALY'], dCOST=m['dCOST'], stategy="Preemptive Panel", fatal_b = 2e+5*m['fatal_b']))

  