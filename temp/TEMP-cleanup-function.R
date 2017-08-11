cleanup_on_termination <- function(traj)
       {
       traj %>% 
       branch(
                                  function(attrs) attrs[["aTreat_SC_A"]]+1,
                                  continue=rep(TRUE,2),
                                  trajectory() %>% timeout(0),
                                  trajectory() %>% 
                                  branch(
                                  function(attrs) attrs[["aDrug_SC_A"]],
                                  continue=rep(TRUE,2),
                                  trajectory() %>% release("rx_SC_A"),
                                  trajectory() %>% release("alt_SC_A")
                                  )) %>% branch(
                                  function(attrs) attrs[["aTreat_SC_B"]]+1,
                                  continue=rep(TRUE,2),
                                  trajectory() %>% timeout(0),
                                  trajectory() %>% 
                                  branch(
                                  function(attrs) attrs[["aDrug_SC_B"]],
                                  continue=rep(TRUE,2),
                                  trajectory() %>% release("rx_SC_B"),
                                  trajectory() %>% release("alt_SC_B")
                                  )) %>% branch(
                                  function(attrs) attrs[["aTreat_SC_C"]]+1,
                                  continue=rep(TRUE,2),
                                  trajectory() %>% timeout(0),
                                  trajectory() %>% 
                                  branch(
                                  function(attrs) attrs[["aDrug_SC_C"]],
                                  continue=rep(TRUE,2),
                                  trajectory() %>% release("rx_SC_C"),
                                  trajectory() %>% release("alt_SC_C")
                                  )) %>% branch(
                                  function(attrs) attrs[["aTreat_SC_D"]]+1,
                                  continue=rep(TRUE,2),
                                  trajectory() %>% timeout(0),
                                  trajectory() %>% 
                                  branch(
                                  function(attrs) attrs[["aDrug_SC_D"]],
                                  continue=rep(TRUE,2),
                                  trajectory() %>% release("rx_SC_D"),
                                  trajectory() %>% release("alt_SC_D")
                                  )) %>%
release("time_in_model")
}