# This file reads in the CSV file of the scenarios.

scenarios <- read.csv("./simple-pgx-scenario-parameters.csv",stringsAsFactors = FALSE) %>% tbl_df(); head(scenarios)
#scenarios <- scenarios %>% select(1:6)
scenario.names <- scenarios %>% select(-param,-type,-value,-range,-description) %>% names()
scenario.ids <- paste0("sc_",letters[1:length(scenario.names)])
scenario.mapping <- cbind.data.frame(scenario.id = scenario.ids,scenarnio.name = scenario.names)
secnario.mapping <- scenario.names 
names(scenario.names) = gsub("sc_","",scenario.ids)

# Get the risk parameters
risks <- scenarios %>% filter(type=="risk") %>% select(-c(2:5))
names(risks) <- c("param",scenario.ids)
risks2 <- risks %>% gather(scenario,value,-param) %>% unite(parameter,c("param","scenario")) %>% data.frame()
risks3 <- as.numeric(risks2$value)
risks.as.list <- setNames(split(risks3, seq(length(risks3))), risks2$parameter)
inputs.main <- append(list(
  vAge = 40,
  vGender = 1,
  vPreemptive = "None", # "None" or "Panel"
  vReactive = "None", # "None" or "Single" or "Panel"
  vHorizon  = 10,
  vN = 100
),risks.as.list)


# Get the disutilities 
disutility <- scenarios %>% filter(type=="disutility") %>% select(-c(2:5))
names(disutility) <- c("param",scenario.ids)
disutility2 <- disutility %>% gather(scenario,value,-param) %>% unite(parameter,c("param","scenario")) %>% data.frame()
disutility3 <- as.numeric(disutility2$value)
disutility.as.list <- setNames(split(disutility3, seq(length(disutility3))), disutility2$parameter)
disutilities = append(list(
  secular_death = 1
),disutility.as.list)

# Get the durations
duration <- scenarios %>% filter(type=="duration") %>% select(-c(2:5))
names(duration) <- c("param",scenario.ids)
duration2 <- duration %>% gather(scenario,value,-param) %>% unite(parameter,c("param","scenario")) %>% data.frame()
duration3 <- as.numeric(duration2$value)
duration.as.list <- setNames(split(duration3, seq(length(duration3))), duration2$parameter)
durations = append(list(
),duration.as.list)

type <- scenarios %>% filter(type=="type") %>% select(-c(2:5))
names(type) <- c("param",scenario.ids)
type2 <- type %>% gather(scenario,value,-param) %>% unite(parameter,c("param","scenario")) %>% data.frame()
type3 <- as.numeric(type2$value)
type.as.list <- setNames(split(type3, seq(length(type3))), type2$parameter)
type = append(list(
  # A = 1,
  # A_c = 0,
  # B_Survive = 0,
  # B_Death = 0,
  secular_death = 0
),type.as.list)

cost <- scenarios %>% filter(type=="cost") %>% select(-c(2:5))
names(cost) <- c("param",scenario.ids)
cost2 <- cost %>% gather(scenario,value,-param) %>% unite(parameter,c("param","scenario")) %>% data.frame()
cost3 <- as.numeric(cost2$value)
cost.as.list <- setNames(split(cost3, seq(length(cost3))), cost2$parameter)
costs = append(list(
  panel_test=250
),cost.as.list)

inputs <- append(inputs.main,list(disutilities=disutilities,durations=durations,type=type,costs=costs))


