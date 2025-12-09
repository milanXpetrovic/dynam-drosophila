library(jsonlite)

convert_to_json <- function(result_element) {
  # Extract relevant fields
  parameters <- result_element$parameters
  standardErrors <- result_element$standardErrors
  logLikelihood <- result_element$logLikelihood
  finalScore <- result_element$finalScore
  finalInformationMatrix <- result_element$finalInformationMatrix
  convergence <- result_element$convergence
  nIterations <- result_element$nIterations
  nEvents <- result_element$nEvents
  formula <- as.character(result_element$formula)  # Convert formula to character
  model <- result_element$model
  subModel <- result_element$subModel
  rightCensored <- result_element$rightCensored
  nParams <- result_element$nParams
  call <- as.character(result_element$call)  # Convert call to character
  coefMat <- result_element$coefMat
  AIC <- result_element$AIC
  BIC <- result_element$BIC

  # Create the JSON-friendly list
  json_data <- list(
    parameters = parameters,
    standardErrors = standardErrors,
    logLikelihood = logLikelihood,
    finalScore = finalScore,
    finalInformationMatrix = finalInformationMatrix,
    convergence = convergence,
    nIterations = nIterations,
    nEvents = nEvents,
    formula = formula,
    model = model,
    subModel = subModel,
    rightCensored = rightCensored,
    nParams = nParams,
    call = call,
    coefMat = coefMat,
    AIC = AIC,
    BIC = BIC
  )

  # Convert the list to JSON string
  json_string <- toJSON(json_data, pretty = TRUE)
  
  return(json_string)
}