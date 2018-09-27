setwd("smdm")

mar <- read.csv("data/markov-icer-cube.csv")
des <- read.csv("data/des-icer-cube.csv")
deq <- read.csv("data/deq-icer-cube.csv")

mar$icer_bias <- deq$ICER - mar$ICER
mar$nmb_bias  <- deq$NMB  - mar$NMB

mar[which(mar$icer_bias == -max(-mar$icer_bias)),]
mar[which(mar$nmb_bias  == -max(-mar$nmb_bias)),]

m <- lm(nmb_bias ~ p_bd + p_g + r_a + r_b + rr_b + c_a + c_bs + c_bd + c_tx + c_alt + c_t + d_a + d_at + d_b + horizon, mar)