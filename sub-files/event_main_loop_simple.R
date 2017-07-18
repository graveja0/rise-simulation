library(simmer)

  ##############################################
 ##
## Helper functions to make code simpler
##
exec <- function(traj, func)
{
  traj %>%
  timeout(function(attrs) {func(attrs); 0})
}

print_attrs <- function(traj)
{
  exec(traj, function(attrs) print(attrs))
}

write_attrs <- function(traj,file=patients)
{
  exec(traj, function(attrs) patients = rbind.fill(file,as.data.frame(attrs)))
}


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
    traj <- set_attribute(traj, event$attr, function(attrs)
    {
      event$time_to_event(attrs, inputs)
    })
  })
  traj
}

# Find the next event based on time
next_event <- function(attrs)
{
  event_time <- Inf
  event      <- NA
  id         <- 0
  for(i in 1:length(event_registry))
  {
    e <- event_registry[[i]]
    tmp_time   <- attrs[[e$attr]]
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
  traj <- timeout(traj, function(attrs)
  {
    # Determine next up
    ne <- next_event(attrs)
    event <-      ne[['event']]
    event_time <- ne[['event_time']]
    
    #cat(" Next up => ",event$name,"\n")
    #cat("            waiting", event_time-now(env), '\n')
    
    # Wait the clock time for the nearest event, minus now()
    event_time - now(env)
  })
  
  # Age them by clock
  traj <- set_attribute(traj,'aAge',function(attrs) attrs[['aAgeInitial']]+(now(env)/365.0))
  
  # Create a handler for every possible event, using their
  # list position as the branch number
  # This will determine the id of the next event
  # Call it's modification function
  # and then update it's next time to event
  args <- lapply(event_registry,FUN=function(e) {
    #print(e$name)   # Good for debugging event loading
    trajectory(e$name) %>%
      #timeout(function(attrs) {cat("executing ",e$name,"\n"); 0}) %>%
      e$func(inputs) %>%
      #timeout(function(attrs) {cat("executed ",e$name,"\n"); 0}) %>%
      set_attribute(e$attr, function(attrs) {now(env)+e$time_to_event(attrs,inputs)})
  })
  args$".trj"    <- traj
  args$option    <- function(attrs) next_event(attrs)$id
  args$continue  <- rep(TRUE,length(event_registry))
  
  traj <- do.call(branch, args)
  
  # Apply reactive events
  lapply(event_registry[sapply(event_registry, function(x) x$reactive)], FUN=function(e){
    traj <- set_attribute(traj, e$attr, function(attrs) {now(env)+e$time_to_event(attrs,inputs)})
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
##
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