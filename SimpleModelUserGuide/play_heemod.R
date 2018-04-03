# library(heemod)
# library(diagram)

###############################
# Question Box:
# Q: one-time cost and temporary disu
# A: add transition states (inconvenient compared to DES)
# Q: time-dependent risk
# A: define_parameters()

# convenience functions to compute transtion probs
# https://pierucci.org/heemod/reference/probability.html
###############################


rm(list=ls())
################################
# simple case: https://pierucci.org/heemod/articles/c_homogeneous.html
# states: A,B,C,D(death)
# strategies: mono,comb
# transition probs: define_transition()
mat_mono <- define_transition(
        .721, .202, .067, .010,
        0,    .581, .407, .012,
        0,    0,    .750, .250,
        0,    0,    0,    1,
        state_names = c("A","B","C","D")
)

# The combined therapy group has its transition probabilities multiplied by rr=0.509
# use the alias C as a convenient way to specify the probability complement, equal to 1−∑Ptrans
rr <- .509
mat_comb <- define_transition(
        C, .202*rr, .067*rr, .010*rr,
        0, C,       .407*rr, .012*rr,
        0, 0,       C,       .250*rr,
        0, 0,       0,       1,
        state_names = c("A","B","C","D")
)

plot(mat_mono)
plot(mat_comb)

# state values:define_state() with cost and effectiveness implications
# both c and e here are per cycle !!!!!
# dispatch_strategy() returns different values depending on the strategy.
cost_zido <- 2278
cost_lami <- 2086
state_A <- define_state(
        cost_health = discount(2756, .06),
        cost_drugs = discount(dispatch_strategy(
                mono = cost_zido,
                comb = cost_zido + cost_lami
        ), .06),
        cost_total = cost_health + cost_drugs,
        life_year = 1
)
state_B <- define_state(
        cost_health = discount(3052, .06),
        cost_drugs = discount(dispatch_strategy(
                mono = cost_zido,
                comb = cost_zido + cost_lami
        ), .06),
        cost_total = cost_health + cost_drugs,
        life_year = 1
)
state_C <- define_state(
        cost_health = discount(9007, .06),
        cost_drugs = discount(dispatch_strategy(
                mono = cost_zido,
                comb = cost_zido + cost_lami
        ), .06),
        cost_total = cost_health + cost_drugs,
        life_year = 1
)
state_D <- define_state(
        cost_health = 0,
        cost_drugs = 0,
        cost_total = 0,
        life_year = 0
)

#strategy definitions: define_strategy() by binding transition probs matrix with defined states
strat_mono <- define_strategy(
        transition = mat_mono,
        state_A,
        state_B,
        state_C,
        state_D
)

strat_comb <- define_strategy(
        transition = mat_comb,
        state_A,
        state_B,
        state_C,
        state_D
)

# run: run_model()
res_mod <- run_model(
        mono = strat_mono,
        comb = strat_comb,
        cycles = 50,
        cost = cost_total,
        effect = life_year
)

#see results
summary(res_mod,
        threshold = c(1000, 5000, 6000, 1e4))



rm(list=ls())
################################
# time-varying case: https://pierucci.org/heemod/articles/d_non_homogeneous.html
# 5 states, 2 strategies
# Death depends on age (our secular death case!!!!!!!)
# Primary THR revision depends on model time (not relevant to our simple case)

# define param for time varying case: define_parameters()
# need to define all params within the call and directly cite the transition probs formula obj in transition matrices
# markov_cycle corresponding to now() in simmer
param <- define_parameters(
        age_init = 60,
        sex = 0,
        
        # age increases with cycles
        age = age_init + markov_cycle,
        
        # operative mortality rates
        omrPTHR = .02, #primary
        omrRTHR = .02, #revision
        
        # re-revision mortality rate
        rrr = .04,
        
        # parameters for calculating primary revision rate
        cons = -5.49094,
        ageC = -.0367,
        maleC = .768536,
        lambda = exp(cons + ageC * age_init + maleC * sex),
        gamma = 1.45367786,
        
        rrNP1 = .260677,
        
        # revision probability of primary procedure
        standardRR = 1 - exp(lambda * ((markov_cycle - 1) ^ gamma -
                                               markov_cycle ^ gamma)),
        np1RR = 1 - exp(lambda * rrNP1 * ((markov_cycle - 1) ^ gamma - 
                                                  markov_cycle ^ gamma)),
        
        # age-related mortality rate (use age-gompertz in our simple case !!!!!!!)
        sex_cat = ifelse(sex == 0, "FMLE", "MLE"),
        mr = get_who_mr(age, sex_cat, country = "GBR", local = TRUE),
        
        # state values
        u_SuccessP = .85,
        u_RevisionTHR = .30,
        u_SuccessR = .75,
        c_RevisionTHR = 5294
)

# transition probs
mat_standard <- define_transition(
        state_names = c(
                "PrimaryTHR",
                "SuccessP",
                "RevisionTHR",
                "SuccessR",
                "Death"
        ),
        0, C, 0,          0, omrPTHR,
        0, C, standardRR, 0, mr,
        0, 0, 0,          C, omrRTHR+mr,
        0, 0, rrr,        C, mr,
        0, 0, 0,          0, 1
)

mat_np1 <- define_transition(
        state_names = c(
                "PrimaryTHR",
                "SuccessP",
                "RevisionTHR",
                "SuccessR",
                "Death"
        ),
        0, C, 0,          0, omrPTHR,
        0, C, np1RR,      0, mr,
        0, 0, 0,          C, omrRTHR+mr,
        0, 0, rrr,        C, mr,
        0, 0, 0,          0, 1
)

# define strategies & states
strat_standard <- define_strategy(
        transition = mat_standard,
        PrimaryTHR = define_state(
                utility = 0,
                cost = 0
        ),
        SuccessP = define_state(
                utility = discount(u_SuccessP, .015),
                cost = 0
        ),
        RevisionTHR = define_state(
                utility = discount(u_RevisionTHR, .015),
                cost = discount(c_RevisionTHR, .06)
        ),
        SuccessR = define_state(
                utility = discount(u_SuccessR, .015),
                cost = 0
        ),
        Death = define_state(
                utility = 0,
                cost = 0
        ),
        starting_values = define_starting_values(
                cost = 394
        )
)

strat_np1 <- define_strategy(
        transition = mat_np1,
        PrimaryTHR = define_state(
                utility = 0,
                cost = 0
        ),
        SuccessP = define_state(
                utility = discount(u_SuccessP, .015),
                cost = 0
        ),
        RevisionTHR = define_state(
                utility = discount(u_RevisionTHR, .015),
                cost = discount(c_RevisionTHR, .06)
        ),
        SuccessR = define_state(
                utility = discount(u_SuccessR, .015),
                cost = 0
        ),
        Death = define_state(
                utility = 0,
                cost = 0
        ),
        starting_values = define_starting_values(
                cost = 579
        )
)

# run the model
res_mod <- run_model(
        standard = strat_standard,
        np1      = strat_np1,
        parameters = param,
        cycles = 60,
        cost = cost,
        effect = utility
)


# incorporate population heterogeneity
# https://pierucci.org/heemod/articles/g_heterogeneity.html

pop <- data.frame(
        age=c(70,80),
        sex=c(1,0),
        .weights=c(80,20)
)

res_h <- update(res_mod, newdata = pop)

res_mod$run_model
res_h$updated_model
res_h$combined_model$run_model