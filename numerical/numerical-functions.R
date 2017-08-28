
draw.parameter.values <- function(scenario.file ="../simple-pgx-scenario-parameters-vogi.csv",PSA.N=10,ws=NULL) {
  # Read in the scenario spreadsheet and map the (long) scenario names to a generic A, B, C, etc.
  if (is.null(ws)) scenarios <- read.csv(scenario.file,stringsAsFactors = FALSE) %>% tbl_df() else
    scenarios <- read_excel(scenario.file,sheet=ws) %>% tbl_df()
  
  scenario.names <- scenarios %>%
    dplyr::select(-param,-type,-value,-psatype,-description,-dplyr::contains("psa_param")) %>% names()
  
  # allow for up to 720 scenarios (just need enough scenario IDs)
  lots.of.letters <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  scenario.ids <- paste0("SC_",lots.of.letters[1:length(scenario.names)])
  scenario.mapping <- cbind.data.frame(scenario.id = scenario.ids,scenarnio.name = scenario.names)
  scenario.mapping <- scenario.names
  names(scenario.names) = gsub("SC_","",scenario.ids)
  
  drawn.parameter.values <- unique(scenarios$type) %>%
    purrr::map(~draw.latin.hypercube(tt=.x,PSA.N=PSA.N,scenarios=scenarios,scenario.names=scenario.names,
                                     scenario.ids=scenario.ids) )
  names(drawn.parameter.values) <- unique(scenarios$type)
  output <- list()
  output[["parameter.draws"]] <- drawn.parameter.values
  output[["scenarios"]] <- scenario.names
  output
}

draw.latin.hypercube <- function(tt,PSA.N,scenarios,scenario.names,scenario.ids) {
  
  require(lhs)
  
  params.full <-
    scenarios %>% filter(type == tt) %>% select(-value, -description)
  names.temp <- names(params.full)
  for (y in scenario.names)
    names.temp <-
    gsub(paste0(paste0("^", y, "_"), "|", paste0("^", y, "$")),
         paste0("SC_", names(scenario.names[which(scenario.names == y)])),
         names.temp)
  names(params.full) <- gsub("psa_param", "_psa", names.temp)
  if (tt != "global")
  {
    params.full %>% reshape2::melt(id.vars = c("param", "psatype", "type")) %>% tbl_df() %>%
      mutate(paramtype = gsub(paste0(scenario.ids, collapse = "|"), "", variable)) %>%
      mutate(
        variable = gsub("_psa1|_psa2", "", variable) ,
        paramtype = gsub("^_", "", paramtype)
      ) %>%
      mutate(paramtype = ifelse(paramtype == "", "value", paramtype)) %>%
      mutate(paramtype = ifelse(
        psatype == "uniform",
        gsub("psa1", "min", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "uniform",
        gsub("psa2", "max", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "beta",
        gsub("psa1", "shape1", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "beta",
        gsub("psa2", "shape2", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "constant",
        gsub("psa1", "constant1", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "constant",
        gsub("psa2", "constant2", paramtype),
        paramtype
      )) %>%
      rename(scenario = variable) %>% reshape2::dcast(param + psatype +
                                                        type + scenario ~ paramtype) %>% data.frame() %>%
      tidyr::unite(parameter, param, scenario) -> params2
  } else
  {
    params.full %>% reshape2::melt(id.vars = c("param", "psatype", "type")) %>%
      mutate(paramtype = gsub(paste0(scenario.ids, collapse ="|"), "", variable)) %>%
      mutate(
        variable = gsub("_psa1|_psa2", "", variable) ,
        paramtype = gsub("^_", "", paramtype)
      ) %>%
      filter(variable == "SC_A") %>%
      mutate(paramtype = ifelse(paramtype == "", "value", paramtype)) %>%
      mutate(paramtype = ifelse(
        psatype == "uniform",
        gsub("psa1", "min", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "uniform",
        gsub("psa2", "max", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "beta",
        gsub("psa1", "shape1", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "beta",
        gsub("psa2", "shape2", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "constant",
        gsub("psa1", "constant1", paramtype),
        paramtype
      )) %>%
      mutate(paramtype = ifelse(
        psatype == "constant",
        gsub("psa2", "constant2", paramtype),
        paramtype
      )) %>%
      rename(scenario = variable) %>% reshape2::dcast(param + type + psatype + scenario ~
                                                        paramtype) %>%
      mutate(scenario = "global")  %>% rename(parameter = param) -> params2
  }
  params.as.list <-
    setNames(split(params2, 1:nrow(params2)), params2$parameter) %>% purrr::map( ~as.list(.x))
  params <- unlist(lapply(params.as.list, function(x)
    x$value))
  
  X <- randomLHS(PSA.N, length(params))
  colnames(X) = names(params)
  
  lhc.draws.transformed <-
    cbind.data.frame(lapply(params.as.list, function(x)
    {
      if (x[["psatype"]] == "beta")
      {
        qbeta(X[, x[["parameter"]]], shape1 = x[["shape1"]], shape2 = x[["shape2"]])
      }
      else if (x[["psatype"]] == "uniform")
      {
        qunif(X[, x[["parameter"]]], min = x[["min"]], max = x[["max"]])
      }
      else if (x[["psatype"]] == "constant")
      {
        rep(x[["value"]], PSA.N)
      }
    }))
  lhc.draws.transformed
}


run.numerical <- function(scenario, ii,parameter.values) {
  d <- model.run(config, ii, "none",times)
  temp <-
    d %>% tbl_df() %>% cbind(as.data.frame(parameter.values$parameter.draws))
  varying1 <-
    suppressWarnings(temp %>% map_dbl(~ var(., na.rm = TRUE)))
  varying <-
    unique(c("dQALY", "dCOST", names(varying1[which(varying1 > 0)])))
  temp <-
    temp[, varying] %>% data.frame() %>% dplyr::mutate(iteration =
                                                         row_number())  %>%
    dplyr::select(iteration,
                  dplyr::contains("QALY"),
                  dplyr::contains("COST"),
                  everything()) %>% mutate(strategy = ss)
  assign(paste0("results.", ss), temp)
  temp
}


add.params <- function(run,parameter.draws,run.name="run") {
  temp <-
    run %>% tbl_df() %>% cbind(as.data.frame(parameter.draws))
  varying1 <-
    suppressWarnings(temp %>% map_dbl(~ var(., na.rm = TRUE)))
  varying <-
    unique(c("dQALY", "dCOST", names(varying1[which(varying1 > 0)])))
  temp <-
    temp[, varying] %>% data.frame() %>% dplyr::mutate(iteration =
                                                         row_number())  %>%
    dplyr::select(iteration,
                  dplyr::contains("QALY"),
                  dplyr::contains("COST"),
                  everything()) %>% mutate(strategy = run.name)
  return(temp)
}
