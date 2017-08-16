risks.as.list <- setNames(split(t(drawn.parameter.values[["risk"]]), seq(nrow(t(drawn.parameter.values[["risk"]])))), colnames(drawn.parameter.values[["risk"]]))
global.as.list <- setNames(split(t(drawn.parameter.values[["global"]]), seq(nrow(t(drawn.parameter.values[["global"]])))), colnames(drawn.parameter.values[["global"]]))

inputs.main <- append(list(
  vAge = 40,
  vGender = 1,
  vPreemptive = "None", # "None" or "Panel"
  vReactive = "None", # "None" or "Single" or "Panel"
  vHorizon  = 10,
  vN = 100
),risks.as.list)
inputs.main <- append(inputs.main,global.as.list)
disutility.as.list <- setNames(split(t(drawn.parameter.values[["disutility"]]), seq(nrow(t(drawn.parameter.values[["disutility"]])))), colnames(drawn.parameter.values[["disutility"]]))
disutilities = append(list(
  secular_death = 1
),disutility.as.list)

duration.as.list <- setNames(split(t(drawn.parameter.values[["duration"]]), seq(nrow(t(drawn.parameter.values[["duration"]])))), colnames(drawn.parameter.values[["duration"]]))
durations = append(list(
),duration.as.list)

type.as.list <- setNames(split(t(drawn.parameter.values[["type"]]), seq(nrow(t(drawn.parameter.values[["type"]])))), colnames(drawn.parameter.values[["type"]]))
type = append(list(
  secular_death = 0
),type.as.list)


cost.as.list <- setNames(split(t(drawn.parameter.values[["cost"]]), seq(nrow(t(drawn.parameter.values[["cost"]])))), colnames(drawn.parameter.values[["cost"]]))
costs = append(list(
  panel_test=drawn.parameter.values$global$panel_test,
  single_test = drawn.parameter.values$global$single_test
),cost.as.list)


inputs <- append(append(inputs.init,inputs.main),list(disutilities=disutilities,durations=durations,type=type,costs=costs))
