library(simmer)


  ##############################################
 ##
## Helper functions for managing counters
##
## Hopefully, no modification required.
##

# Create the counters, takes a list
create_counters <- function(env, counters)
{
  sapply(counters, FUN=function(counter)
  {
    env <- add_resource(env, counter, Inf, 0)
  })
  
  env
}

# Mark a counter
mark <- function(traj, counter)
{
  traj               %>%
    seize(counter,1)   %>%
    timeout(0)         %>%
    release(counter,1)
}

  ##############################################
 ##
## Helper functions for managing events
##
## Hopefully, no modification required.
##
assign_events <- function(traj, inputs)
{
  sapply(event_registry, FUN=function(event)
  {
    traj <- set_attribute(traj, event$attr, function()
    {
      event$time_to_event(inputs)
    })
  })
  traj
}

# Find the next event based on time
next_event <- function()
{
  event_time <- Inf
  event      <- NA
  id         <- 0
  for(i in 1:length(event_registry))
  {
    e <- event_registry[[i]]
    tmp_time   <- get_attribute(env,e$attr)
    if(tmp_time < event_time)
    {
      event      <- e
      event_time <- tmp_time
      id         <- i
    }
  }
  
  return(list(event=event, event_time=event_time, id=id))
}

# Process events in main loop
process_events <- function(traj, env, inputs)
{
  # Find the next event from possible events, and timeout (wait) till that moment
  traj <- timeout(traj, function()
  {
    # Determine next up
    ne <- next_event()
    event      <- ne[['event']]
    event_time <- ne[['event_time']]
    
    #cat(" Next up => ",event$name,"\n")
    #cat("            waiting", event_time-now(env), '\n')
    
    # Wait the clock time for the nearest event, minus now()
    event_time - now(env)
  })
  
  # Age them by clock
  traj <- set_attribute(traj,'aAge',function() get_attribute(env,"aAgeInitial")+(now(env)/365.0))
  
  # Create a handler for every possible event, using their
  # list position as the branch number
  # This will determine the id of the next event
  # Call it's modification function
  # and then update it's next time to event
  args <- lapply(event_registry,FUN=function(e) {
    #print(e$name)   # Good for debugging event loading
    trajectory(e$name) %>%
      e$func(inputs) %>%
      set_attribute(e$attr, function() {now(env)+e$time_to_event(inputs)})
  })
  args$".trj"    <- traj
  args$option    <- function() next_event()$id
  args$continue  <- rep(TRUE,length(event_registry))
  
  traj <- do.call(branch, args)
  
  # Apply reactive events
  lapply(event_registry[sapply(event_registry, function(x) x$reactive)], FUN=function(e){
    traj <- set_attribute(traj, e$attr, function() {now(env)+e$time_to_event(inputs)})
  })
  
  traj
}

  ##############################################
 ##
## MAIN LOOP
##
## This should not require modification
## This creates a patient simulation (trajectory)
## 
## It uses a branch in a manner to prevent the
## rollback from looking further up the stack
## of the event loop. 
simulation <- function(env, inputs)
{
  trajectory("Patient")     %>%
    initialize_patient(inputs)     %>%
    assign_events(inputs)          %>%
    preemptive_strategy(inputs)    %>% 
    branch( # Used branch, to prevent rollback from looking inside event loop function
      function() 1,
      continue=TRUE,
      trajectory("main_loop") %>% process_events(env, inputs)
    ) %>% 
    rollback(amount=1, times=100) # Process up to 100 events per person
}