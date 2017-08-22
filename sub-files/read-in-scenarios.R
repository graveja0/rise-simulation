

# Read in the scenario spreadsheet and map the (long) scenario names to a generic A, B, C, etc.
scenarios <- read.csv(scenario.file,stringsAsFactors = FALSE) %>% tbl_df(); head(scenarios)

#run.psa <- FALSE
#if (!(run.psa)) scenarios$psatype = "constant"
scenario.names <- scenarios %>% dplyr::select(-param,-type,-value,-psatype,-description,-dplyr::contains("psa_param")) %>% names()

# allow for up to 720 scenarios (just need enough scenario IDs)
lots.of.letters <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
scenario.ids <- paste0("SC_",lots.of.letters[1:length(scenario.names)])
scenario.mapping <- cbind.data.frame(scenario.id = scenario.ids,scenarnio.name = scenario.names)
scenario.mapping <- scenario.names
names(scenario.names) = gsub("SC_","",scenario.ids)

tt <- "risk"

require(lhs)
require(tidyverse)

draw.latin.hypercube <- function(tt,PSA.N=10) {
  params.full <- scenarios %>% filter(type==tt) %>% select(-value,-description)
  names.temp <- names(params.full)
  #for (y in scenario.names) names.temp <- gsub(paste0("^",y,"_"),paste0("SC_",names(scenario.names[which(scenario.names==y)])),names.temp)
  for (y in scenario.names) names.temp <- gsub( paste0(paste0("^",y,"_"),"|",paste0("^",y,"$")),paste0("SC_",names(scenario.names[which(scenario.names==y)])),names.temp)
  names(params.full) <- gsub("psa_param","_psa",names.temp)
  if (tt!="global")
  {
    params.full %>% melt(id.vars = c("param","psatype","type")) %>% tbl_df() %>% 
      mutate(paramtype = gsub(paste0(scenario.ids,collapse="|"),"",variable)) %>%
      mutate(variable = gsub("_psa1|_psa2","",variable) , paramtype = gsub("^_","",paramtype)) %>%
      mutate(paramtype = ifelse(paramtype=="","value",paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="uniform",gsub("psa1","min",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="uniform",gsub("psa2","max",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="beta",gsub("psa1","shape1",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="beta",gsub("psa2","shape2",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="constant",gsub("psa1","constant1",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="constant",gsub("psa2","constant2",paramtype),paramtype)) %>%
      rename(scenario = variable) %>% reshape2::dcast(param+psatype+type+scenario~paramtype) %>% data.frame() %>% 
      tidyr::unite(parameter,param,scenario) -> params2
  } else 
  {
    params.full %>% melt(id.vars = c("param","psatype","type")) %>% mutate(paramtype = gsub(paste0(scenario.ids,collapse="|"),"",variable)) %>%
      mutate(variable = gsub("_psa1|_psa2","",variable) , paramtype = gsub("^_","",paramtype)) %>%
      filter(variable=="SC_A") %>% 
      mutate(paramtype = ifelse(paramtype=="","value",paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="uniform",gsub("psa1","min",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="uniform",gsub("psa2","max",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="beta",gsub("psa1","shape1",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="beta",gsub("psa2","shape2",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="constant",gsub("psa1","constant1",paramtype),paramtype)) %>%
      mutate(paramtype = ifelse(psatype=="constant",gsub("psa2","constant2",paramtype),paramtype)) %>%
      rename(scenario = variable) %>% dcast(param+type+psatype+scenario~paramtype) %>% 
      mutate(scenario = "global")  %>% rename(parameter = param)-> params2
  }
  params.as.list <- setNames(split(params2, 1:nrow(params2)), params2$parameter) %>% purrr::map(~as.list(.x))
  params <- unlist(lapply(params.as.list,function(x) x$value))

  X <- randomLHS(inputs.init$vN_PSA, length(params))
  colnames(X) = names(params)

  lhc.draws.transformed <- cbind.data.frame(lapply(params.as.list,function(x)
  {
    if (x[["psatype"]]=="beta")
    {
      qbeta(X[,x[["parameter"]]],shape1=x[["shape1"]],shape2=x[["shape2"]])
    }
    else if (x[["psatype"]]=="uniform")
    {
      qunif(X[,x[["parameter"]]],min=x[["min"]],max=x[["max"]])
    }
    else if (x[["psatype"]]=="constant")
    {
      rep(x[["value"]],PSA.N)
    }
  }
  ))
  lhc.draws.transformed
}
#ii = 1

drawn.parameter.values <- unique(scenarios$type) %>% purrr::map(~draw.latin.hypercube(tt=.x,PSA.N=inputs.init$vN_PSA) )
names(drawn.parameter.values) <- unique(scenarios$type)


# risks.as.list <- setNames(split(t(drawn.parameter.values[["risk"]][ii,]), seq(nrow(t(drawn.parameter.values[["risk"]][ii,])))), colnames(drawn.parameter.values[["risk"]]))
# inputs.main <- append(list(
#   vAge = 40,
#   vGender = 1,
#   vPreemptive = "None", # "None" or "Panel"
#   vReactive = "None", # "None" or "Single" or "Panel"
#   vHorizon  = 10,
#   vN = 100
# ),risks.as.list)
# disutility.as.list <- setNames(split(t(drawn.parameter.values[["disutility"]][ii,]), seq(nrow(t(drawn.parameter.values[["disutility"]][ii,])))), colnames(drawn.parameter.values[["disutility"]]))
# disutilities = append(list(
#   secular_death = 1
# ),disutility.as.list)
# 
# duration.as.list <- setNames(split(t(drawn.parameter.values[["duration"]][ii,]), seq(nrow(t(drawn.parameter.values[["duration"]][ii,])))), colnames(drawn.parameter.values[["duration"]]))
# durations = append(list(
# ),duration.as.list)
# 
# type.as.list <- setNames(split(t(drawn.parameter.values[["type"]][ii,]), seq(nrow(t(drawn.parameter.values[["type"]][ii,])))), colnames(drawn.parameter.values[["type"]]))
# type = append(list(
#   secular_death = 0
# ),type.as.list)
# 
# 
# cost.as.list <- setNames(split(t(drawn.parameter.values[["cost"]][ii,]), seq(nrow(t(drawn.parameter.values[["cost"]][ii,])))), colnames(drawn.parameter.values[["cost"]]))
# costs = append(list(
#   panel_test=250
# ),cost.as.list)
# 
# inputs <- append(inputs.main,list(disutilities=disutilities,durations=durations,type=type,costs=costs))

