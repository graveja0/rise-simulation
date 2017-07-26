risks.as.list <- setNames(split(t(drawn.parameter.values[["risk"]]), seq(nrow(t(drawn.parameter.values[["risk"]])))), colnames(drawn.parameter.values[["risk"]]))
inputs.main <- append(list(
  vAge = 40,
  vGender = 1,
  vPreemptive = "None", # "None" or "Panel"
  vReactive = "None", # "None" or "Single" or "Panel"
  vHorizon  = 10,
  vN = 100
),risks.as.list)
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
  panel_test=250
),cost.as.list)

inputs <- append(append(inputs.init,inputs.main),list(disutilities=disutilities,durations=durations,type=type,costs=costs))
