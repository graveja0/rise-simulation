## Secular Death, Weibull Model
source('./age-gompertz.R')

# Given attributes of a patient (trajectory), it returns in days 
# how long till the patient would die a secular death.
#
# NOTE: The variable in must be named attrs
days_till_death <- function(attrs, inputs)
{
  age       <- attrs[['aAge']]
  death_age <- ageAtDeath(age, attrs[['aGender']])
  
  return(365*(death_age-age))
}

# Given a trajectory, modify as needed when a secular
# death occurs.
#
# In this case, it marks a counter and terminates 
# the trajectory. A branch is required, even though
# it doesn't branch to force the termination.
secular_death <- function(traj, inputs)
{
  traj %>% branch(
    function() 1,
    continue=c(FALSE), # False is patient death, had to use a branch to force termination
    trajectory("Secular Death") %>% mark("secular_death") %>% terminate_simulation()
  )
}