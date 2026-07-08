# ------------------------------------------------------------------
# This file contains helper functions for the DyNAM analysis.
# To use these functions in another script, use: source("utils.R")
# ------------------------------------------------------------------

# Required libraries for the functions in this file
library(jsonlite)
library(dplyr)
library(stringr)

# --- Function Definitions ---

#' Gathers data from multiple goldfish group objects into a single list
#' suitable for stan.
#'
#' @param groupsInfo A data frame with information about each group.
#' @param data A list of data objects, one for each group.
#' @param keepChoice An optional character vector of effect names to keep.
#' @return A gathered data object.
GatherGroups <- function(
    groupsInfo,
    data,
    keepChoice = NULL
) {
  stopifnot(
    is.list(data),
    is.data.frame(groupsInfo),
    is.null(keepChoice) || is.character(keepChoice),
    length(data) == nrow(groupsInfo)
  )
  # init data objects
  colNamesX <- colnames(data[[1]]$data_stan$X_choice)
  
  Nchoice <- sum(sapply(data, \(x) x$data_stan$N_choice))
  Tchoice <- sum(sapply(data, \(x) x$data_stan$T_choice))
  Pchoice <- max(sapply(data, \(x) x$data_stan$P_choice))
  if (!is.null(keepChoice)) {
    if (!all(keepChoice %in% colNamesX))
      stop("some vars in the argument 'keepChoice' are not in the matrix")
    
    Pchoice <- length(keepChoice)
  }
  Xchoice <- matrix(
    0, nrow = Nchoice, ncol = Pchoice,
    dimnames = list(seq_len(Nchoice), colNamesX)
  )
  Qchoice <- max(sapply(data, \(x) x$data_stan$Q_choice))
  # Zchoice <- if (Qchoice == 1) numeric(Nchoice) else
  #   matrix(0, nrow = Nchoice, ncol = Qchoice,
  #          dimnames = list(seq_len(Nchoice),
  #                          colnames(data[[1]]$data_stan$Z_choice)))
  Zchoice <- matrix(0, nrow = Nchoice, ncol = Qchoice,
                  dimnames = list(NULL, colnames(data[[1]]$data_stan$Z_choice)))
  
  nA <- nrow(groupsInfo)
  
  startChoice <- integer(Tchoice)
  endChoice <- integer(Tchoice)
  choseChoice <- integer(Tchoice)
  senderChoice <- integer(Nchoice)
  start_group <- integer(nA)

  if (Pchoice == 1)
    warning("some errors might happen due to only use one covariate")
  
  N <- 1
  E <- 1
  
  for (gr in seq_len(nrow(groupsInfo))) {
    Ngr <- N + data[[gr]]$data_stan$N_choice - 1
    Egr <- E + data[[gr]]$data_stan$T_choice - 1

    Xaux <- data[[gr]]$data_stan$X_choice
    currentColnames <- colnames(Xaux)
    if (Pchoice != length(currentColnames))
      stop("dimension of X doesn't conform for group: ", groupsInfo$group[gr],
           " previous groups ncol: ", P, " this group: ",
           length(currentColnames))
    if (!identical(colNamesX, currentColnames))
      warning("colum names from previous groups is not the same as current",
              " for group: ", groupsInfo$group[gr], "\n\t previous names: ",
              paste(colNamesX, collapse = ", "), "\n\t current names: ",
              paste(currentColnames, collapse = ", "))
    
    if (!is.null(keepChoice)) {
      Xaux <- Xaux[, keepChoice]
    }
    
    posX <- seq(N, Ngr)
    
    Xchoice[posX, ] <- Xaux
    
    if (Qchoice == 1) {
      Zchoice[posX] <- data[[gr]]$data_stan$Z_choice
    } else Zchoice[posX, ] <- data[[gr]]$data_stan$Z_choice
    
    posE <- seq(E, Egr)
    startChoice[posE] <- data[[gr]]$data_stan$start_choice + N - 1
    endChoice[posE] <- data[[gr]]$data_stan$end_choice + N - 1
    choseChoice[posE] <- data[[gr]]$data_stan$chose_choice + N - 1
    senderChoice[posX] <- groupsInfo$ixGroup[gr]
    
    # update N and E
    N <- Ngr + 1
    E <- Egr + 1

    start_group[gr] <- posE[1]

  }
  
  dataGathered <- list(
    data_stan = list(
      T_choice = Tchoice,
      N_choice = Nchoice,
      P_choice = Pchoice,
      Q_choice = data[[1]]$data_stan$Q_choice,
      Qchoice = data[[1]]$data_stan$Q_choice, ##! left this for make_model_code
      A = nA,
      start_choice = startChoice,
      end_choice = endChoice,
      X_choice = Xchoice,
      Z_choice = Zchoice,
      chose_choice = choseChoice,
      sender = senderChoice,
      start_group = start_group
    ),
    namesEffects = data[[1]]$namesEffects,
    effectsDescription = data[[1]]$effectsDescription,
    groupInfo = groupsInfo
  )
  class(dataGathered) <- class(data[[1]])
  attr(dataGathered, "model") <- attr(data[[1]], "model")
  attr(dataGathered, "subModel") <- attr(data[[1]], "subModel")
  return(dataGathered)
}


#' Extracts a unique, sorted list of actors from an event data frame.
#'
#' @param df A data frame with 'sender' and 'receiver' columns.
#' @return A data frame with a single 'label' column of actor names.
getActors <- function(df) {
  unique_actors <- sort(unique(c(df$sender, df$receiver)))
  unique_actors <- unique_actors[unique_actors != "all"]
  unique_actors_df <- data.frame(label = unique_actors)
  return(unique_actors_df)
}


#' Converts a goldfish summary object to a JSON string.
#'
#' @param result_element The summary object from a goldfish model.
#' @return A JSON formatted string.
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

  json_string <- toJSON(json_data, pretty = TRUE)
  return(json_string)
}

#' Reads and cleans a covariate CSV file.
#'
#' @param filename The name of the CSV file.
#' @param group_path The path to the directory containing the file.
#' @return A cleaned data frame.
read_and_clean <- function(filename, group_path) {
  # Note: I modified this function slightly to be more general
  # by accepting group_path as an argument.
  data <- read.csv(file.path(group_path, filename))
  if (ncol(data) == 4) data <- data[, -1]
  colnames(data) <- c("time", "node", "replace")
  return(data)
}

### refractored foo
gather_group_data_for_stan <- function(
    group_metadata,
    group_data_list,
    effects_to_keep = NULL
) {
  stopifnot(
    is.list(group_data_list),
    is.data.frame(group_metadata),
    is.null(effects_to_keep) || is.character(effects_to_keep),
    length(group_data_list) == nrow(group_metadata)
  )

  first_group_data <- group_data_list[[1]]$data_stan
  model_effect_names <- colnames(first_group_data$X_choice)

  total_events <- sum(sapply(group_data_list, \(d) d$data_stan$T_choice))
  total_choice_options <- sum(sapply(group_data_list, \(d) d$data_stan$N_choice))

  num_x_predictors <- max(sapply(group_data_list, \(d) d$data_stan$P_choice))
  num_z_predictors <- max(sapply(group_data_list, \(d) d$data_stan$Q_choice))

  if (!is.null(effects_to_keep)) {
    if (!all(effects_to_keep %in% model_effect_names)) {
      stop("Some variable names in 'effects_to_keep' are not in the design matrix.")
    }
    num_x_predictors <- length(effects_to_keep)
  }

  num_groups <- nrow(group_metadata)

  final_x_colnames <- if (!is.null(effects_to_keep)) effects_to_keep else model_effect_names
  
  combined_X_matrix <- matrix(
    0,
    nrow = total_choice_options,
    ncol = num_x_predictors,
    dimnames = list(NULL, final_x_colnames)
  )

  combined_Z_matrix <- matrix(
    0,
    nrow = total_choice_options,
    ncol = num_z_predictors,
    dimnames = list(NULL, colnames(first_group_data$Z_choice))
  )

  event_start_indices <- integer(total_events)
  event_end_indices <- integer(total_events)
  chosen_option_indices <- integer(total_events)
  option_group_id <- integer(total_choice_options)
  group_start_event_index <- integer(num_groups)

  if (num_x_predictors == 1) {
    warning("Using only one covariate may cause issues with some models or matrix operations.")
  }

  choice_row_cursor <- 1
  event_cursor <- 1

  for (group_idx in seq_len(num_groups)) {
    current_group_stan_data <- group_data_list[[group_idx]]$data_stan

    group_end_choice_row <- choice_row_cursor + current_group_stan_data$N_choice - 1
    group_end_event <- event_cursor + current_group_stan_data$T_choice - 1

    current_choice_rows <- seq(choice_row_cursor, group_end_choice_row)
    current_event_range <- seq(event_cursor, group_end_event)

    current_group_X <- current_group_stan_data$X_choice
    current_colnames <- colnames(current_group_X)

    if (ncol(current_group_X) != length(model_effect_names)) {
      stop(
        "Dimension of X matrix for group: ", group_metadata$group[group_idx],
        " does not conform. Expected ", length(model_effect_names),
        " columns, but got ", ncol(current_group_X), "."
      )
    }
    if (!identical(model_effect_names, current_colnames)) {
      warning(
        "Column names for group: ", group_metadata$group[group_idx],
        " differ from the first group."
      )
    }

    if (!is.null(effects_to_keep)) {
      current_group_X <- current_group_X[, effects_to_keep, drop = FALSE]
    }

    combined_X_matrix[current_choice_rows, ] <- current_group_X
    combined_Z_matrix[current_choice_rows, ] <- current_group_stan_data$Z_choice

    offset <- choice_row_cursor - 1
    event_start_indices[current_event_range] <- current_group_stan_data$start_choice + offset
    event_end_indices[current_event_range] <- current_group_stan_data$end_choice + offset
    chosen_option_indices[current_event_range] <- current_group_stan_data$chose_choice + offset

    option_group_id[current_choice_rows] <- group_metadata$ixGroup[group_idx]
    group_start_event_index[group_idx] <- current_event_range[1]

    choice_row_cursor <- group_end_choice_row + 1
    event_cursor <- group_end_event + 1
  }

  gathered_data <- list(
    data_stan = list(
      T_choice = total_events,
      N_choice = total_choice_options,
      P_choice = num_x_predictors,
      Q_choice = num_z_predictors,
      A = num_groups,
      start_choice = event_start_indices,
      end_choice = event_end_indices,
      X_choice = combined_X_matrix,
      Z_choice = combined_Z_matrix,
      chose_choice = chosen_option_indices,
      sender = option_group_id,
      start_group = group_start_event_index
    ),
    namesEffects = group_data_list[[1]]$namesEffects,
    effectsDescription = group_data_list[[1]]$effectsDescription,
    groupInfo = group_metadata
  )

  class(gathered_data) <- class(group_data_list[[1]])
  attr(gathered_data, "model") <- attr(group_data_list[[1]], "model")
  attr(gathered_data, "subModel") <- attr(group_data_list[[1]], "subModel")

  return(gathered_data)
}


gather_group_data_for_stan_vectorized <- function(
    group_metadata,
    group_data_list,
    effects_to_keep = NULL
) {
  stopifnot(
    is.list(group_data_list),
    is.data.frame(group_metadata),
    is.null(effects_to_keep) || is.character(effects_to_keep),
    length(group_data_list) == nrow(group_metadata)
  )

  stan_data_list <- lapply(group_data_list, `[[`, "data_stan")
  
  n_choices_per_group <- vapply(stan_data_list, `[[`, "N_choice", FUN.VALUE = integer(1))
  n_events_per_group <- vapply(stan_data_list, `[[`, "T_choice", FUN.VALUE = integer(1))
  
  first_group_data <- stan_data_list[[1]]
  model_effect_names <- colnames(first_group_data$X_choice)
  
  combined_X_matrix <- do.call(rbind, lapply(stan_data_list, `[[`, "X_choice"))
  combined_Z_matrix <- do.call(rbind, lapply(stan_data_list, `[[`, "Z_choice"))
  
  if (!is.null(effects_to_keep)) {
    if (!all(effects_to_keep %in% model_effect_names)) {
      stop("Some variable names in 'effects_to_keep' are not in the design matrix.")
    }
    combined_X_matrix <- combined_X_matrix[, effects_to_keep, drop = FALSE]
  }

  offsets <- cumsum(c(0, n_choices_per_group[-length(n_choices_per_group)]))
  
  add_offset <- function(vec, offset) vec + offset
  event_start_indices <- unlist(Map(add_offset, lapply(stan_data_list, `[[`, "start_choice"), offsets))
  event_end_indices   <- unlist(Map(add_offset, lapply(stan_data_list, `[[`, "end_choice"),   offsets))
  chosen_option_indices <- unlist(Map(add_offset, lapply(stan_data_list, `[[`, "chose_choice"), offsets))
  
  option_group_id <- rep(group_metadata$ixGroup, times = n_choices_per_group)
  group_start_event_index <- cumsum(c(1, n_events_per_group[-length(n_events_per_group)]))

  gathered_data <- list(
    data_stan = list(
      T_choice = sum(n_events_per_group),
      N_choice = sum(n_choices_per_group),
      P_choice = ncol(combined_X_matrix),
      Q_choice = ncol(combined_Z_matrix),
      A = nrow(group_metadata),
      start_choice = event_start_indices,
      end_choice = event_end_indices,
      X_choice = combined_X_matrix,
      Z_choice = combined_Z_matrix,
      chose_choice = chosen_option_indices,
      sender = option_group_id,
      start_group = group_start_event_index
    ),
    namesEffects = group_data_list[[1]]$namesEffects,
    effectsDescription = group_data_list[[1]]$effectsDescription,
    groupInfo = group_metadata
  )
  
  class(gathered_data) <- class(group_data_list[[1]])
  attr(gathered_data, "model") <- attr(group_data_list[[1]], "model")
  attr(gathered_data, "subModel") <- attr(group_data_list[[1]], "subModel")
  
  return(gathered_data)
}