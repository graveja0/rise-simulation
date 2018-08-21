rm(list=ls())

x <- read.csv("psid_ch13.csv")

N <- 100


inh_rate <- 0.5 # rate of 1 inheritance 1 family unit away

# Let's assume everyone in family set was diagnosed with some genetic condition and
# estimate propagation rates

sim <- do.call(rbind, lapply(1:N, function(i) {
  y <- x[i,]
  result <- data.frame()
  if(!is.na(y$agem) && y$agem > 0)
    result <- rbind(result, data.frame(sex="F", age=y$agem, prob_d=inh_rate))
  if(!is.na(y$agef) && y$agef > 0)
    result <- rbind(result, data.frame(sex="M", age=y$agef, prob_d=inh_rate))
  if(!is.na(y$num_chd) && y$num_chd > 0)
  {
    for(j in 1:y$num_chd)
    {
      age <- y[,paste0("age_chd",j)]
      sex <- y[,paste0("sex_chd",j)]
      if(!is.na(sex) && !is.na(age))
      {
        result <- rbind(result, data.frame(sex=ifelse(sex == 1, "M", "F"), age=age, prob_d=inh_rate))
      }
    }
  }
  
  if(!is.na(y$num_sib) && y$num_sib > 0)
  {
    for(j in 1:y$num_sib)
    {
      age <- y[,paste0("age_sib",j)]
      sex <- y[,paste0("sex_sib",j)]
      typ <- y[,paste0("sibtype",j)]
      if(!is.na(sex) && !is.na(age) && !is.na(typ) && typ > 0 && typ < 6)
      {
        rate <- if(typ == 1)     inh_rate else
                if(typ == 2) 0.5*inh_rate else
                if(typ == 3) 0.5*inh_rate else
                if(typ == 4) 0.5*inh_rate else
                             0.5*inh_rate
        result <- rbind(result, data.frame(sex=ifelse(sex == 1, "M", "F"), age=age, prob_d=rate))
      }
    }
  }
  
  # Finish up with source
  if(length(result$sex) > 0)
  {
    result$src_row <- i
    result$src_sex <- ifelse(y$sex == 1, "M", "F")
    result$src_age <- y$agei
  }
  
  result
}))






## This was to reproduce structure of family data
## Pursuing estimation via straight sampling
# N <- 1000
# 
# m_ref <- subset(x, sex==1)
# f_ref <- subset(x, sex==2)
# 
# # Cleanup of age data (i.e. 0 is not a valid age for a parent)
# # And the child has to be older than the parent
# m_ref$agef <- ifelse(m_ref$agef <= 0 | (m_ref$agef - 12) < m_ref$agei, NA, m_ref$agef)
# m_ref$agem <- ifelse(m_ref$agem <= 0 | (m_ref$agem - 12) < m_ref$agei, NA, m_ref$agem)
# f_ref$agef <- ifelse(f_ref$agef <= 0 | (f_ref$agef - 12) < f_ref$agei, NA, f_ref$agef)
# f_ref$agem <- ifelse(f_ref$agem <= 0 | (f_ref$agem - 12) < f_ref$agei, NA, f_ref$agem)
# 
# # Now first step is to determine sex based on probability in original population
# m_prob <- dim(m_ref)[1] / (dim(m_ref)[1] + dim(f_ref)[1])
# sim <- data.frame(gender = sample(factor(c("M","F")), N, replace=TRUE, prob=c(m_prob, 1-m_prob)))
# 
# # Next, based on sex pull individual age from original reference distributions
# sim$agei <- ifelse(sim$gender == "M",
#                    sample(m_ref$agei[!is.na(m_ref$agei)], N, replace=TRUE),
#                    sample(f_ref$agei[!is.na(f_ref$agei)], N, replace=TRUE))
# 
# # Determine parent's age from regression model with noise
# m_reg <- lm(agem ~ agei, m_ref)
# f_reg <- lm(agem ~ agei, f_ref)
# sim$agem <- ifelse(sim$gender == "M",
#                    coef(m_reg)[1] + coef(m_reg)[2]*sim$agei+rnorm(N, sd=sd(m_reg$residuals)),
#                    coef(f_reg)[1] + coef(f_reg)[2]*sim$agei+rnorm(N, sd=sd(f_reg$residuals)))
# m_reg <- lm(agef ~ agei, m_ref)
# f_reg <- lm(agef ~ agei, f_ref)
# sim$agef <- ifelse(sim$gender == "M",
#                    coef(m_reg)[1] + coef(m_reg)[2]*sim$agei+rnorm(N, sd=sd(m_reg$residuals)),
#                    coef(f_reg)[1] + coef(f_reg)[2]*sim$agei+rnorm(N, sd=sd(f_reg$residuals)))
# 
# plot(m_ref$agei, m_ref$agef, xlab="Age", ylab="Father's Age")
# points(sim$agei, sim$agef, col='red', pch=2)
# 
# # Number of children from negative binomial regression
# # Ref: http://data.library.virginia.edu/getting-started-with-negative-binomial-regression-modeling/
# library(MASS)
# m_reg <- glm.nb(num_chd ~ agei, m_ref)
# f_reg <- glm.nb(num_chd ~ agei, f_ref)
# sim$num_chd <- ifelse(sim$gender=="M",
#                       rnbinom(N, size=m_reg$theta, mu = exp(sum(coef(m_reg)))),
#                       rnbinom(N, size=f_reg$theta, mu = exp(sum(coef(f_reg)))))
  
                      