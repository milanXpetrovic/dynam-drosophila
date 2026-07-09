#!/usr/bin/env Rscript
# =============================================================================
# RATE stage 1 of 2 — PREPROCESS.  Run ONCE:
#   cd /home/milky/droso-pipe/6_dynam_models
#   nohup Rscript rate_1_prep.r > res/rate/prep.log 2>&1 &
#   tail -f res/rate/prep.log
# Produces res/rate/groups_data_rate.rds + groups_info_rate.rds.
# Then run rate_2_fit.r (which reloads these from disk).
# =============================================================================

suppressPackageStartupMessages({
  library(goldfish.latent); library(goldfish); library(glue)
  library(foreach); library(doParallel); library(parallel)
})

# ===== SET THESE ONCE =====
PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
SRC_DIR <- file.path(PROJECT, "src")

base_path        <- PROJECT
edgelist_path    <- file.path(base_path, "data", "edgelists")
covariances_path <- file.path(base_path, "data", "covariances")
out_path         <- file.path(base_path, "res", "rate")
utils_file       <- file.path(SRC_DIR, "utils.R")
treatments_list  <- c("CS_10D", "Cs_5DIZ", "CsCh")

if (!dir.exists(edgelist_path)) stop("No data/edgelists -> fix PROJECT.\n  ", edgelist_path)
if (!file.exists(utils_file))   stop("utils.R not found -> fix SRC_DIR.\n  ", utils_file)
source(utils_file)
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

run_log   <- file.path(out_path, format(Sys.time(), "prep_%Y%m%d_%H%M%S.log"))
group_log <- file.path(out_path, "step1_groups.log"); file.create(group_log, showWarnings = FALSE)
log_msg <- function(...) { line <- paste0(format(Sys.time(), "%F %T"), "  ", ...)
  cat(line, "\n", sep = "", file = run_log, append = TRUE); message(line) }
log_msg("=== rate prep START ===")

rate_model_formula <- dependent ~ 1 +
  indeg(interaction_network, weighted = FALSE, window = 96) +
  indeg(interaction_network, weighted = FALSE, window = 288) +
  indeg(interaction_network, weighted = FALSE, window = 864) +
  indeg(interaction_network, weighted = FALSE, window = 2592) +
  outdeg(interaction_network, weighted = FALSE, window = 96) +
  outdeg(interaction_network, weighted = FALSE, window = 288) +
  outdeg(interaction_network, weighted = FALSE, window = 864) +
  outdeg(interaction_network, weighted = FALSE, window = 2592) +
  ego(nodesAttr$positive_influence) +
  ego(nodesAttr$negative_influence) +
  ego(nodesAttr$positive_inf_weighted) +
  ego(nodesAttr$negative_inf_weighted) +
  ego(nodesAttr$distance_traveled_between_interactions) +
  ego(nodesAttr$number_of_flies_in_soc_space) +
  ego(nodesAttr$unique_partners_met_social_space) +
  ego(flies$young) + ego(flies$old) + ego(flies$isolated)

all_files_list <- lapply(treatments_list, function(treatment) {
  csv_files <- list.files(file.path(edgelist_path, treatment), pattern = "*.csv", full.names = TRUE)
  if (length(csv_files) > 0)
    data.frame(treatment_name = treatment, file_path = csv_files,
               group_name = tools::file_path_sans_ext(basename(csv_files)))
  else NULL
})
files_to_process <- do.call(rbind, all_files_list)
log_msg("groups found: ", nrow(files_to_process))

cl <- makeCluster(parallel::detectCores() - 2)
registerDoParallel(cl)

parallel_results <- foreach(
  i = seq_len(nrow(files_to_process)),
  .packages = c("goldfish", "goldfish.latent", "glue", "tools"),
  .errorhandling = "pass"
) %dopar% {
  treatment_name <- files_to_process$treatment_name[i]
  group_name     <- files_to_process$group_name[i]
  .t0 <- Sys.time()
  cat(sprintf("%s [PID %d] START %2d/%d  %s\n", format(.t0, "%F %T"), Sys.getpid(),
              i, nrow(files_to_process), group_name), file = group_log, append = TRUE)

  interaction_data <- read.csv(files_to_process$file_path[i])
  interaction_data <- interaction_data[, c("time", "sender", "receiver", "increment")]

  actors <- sort(unique(c(interaction_data$sender, interaction_data$receiver)))
  flies  <- data.frame(label = as.character(actors), present = TRUE, stringsAsFactors = FALSE)
  for (col in c("popularity", "activity", "popularity_weighted", "activity_weighted",
                "positive_influence", "negative_influence",
                "positive_inf_weighted", "negative_inf_weighted",
                "distance_traveled_between_interactions", "number_of_flies_in_soc_space",
                "unique_partners_met_interaction_space", "unique_partners_met_social_space",
                "young", "old", "isolated")) flies[[col]] <- 0
  if (treatment_name == "CsCh")    flies$young    <- 1
  if (treatment_name == "CS_10D")  flies$old      <- 1
  if (treatment_name == "Cs_5DIZ") flies$isolated <- 1

  gcp <- file.path(covariances_path, treatment_name, group_name)
  cov_in_degree            <- read_and_clean("in_degree.csv", gcp)
  cov_out_degree           <- read_and_clean("out_degree.csv", gcp)
  cov_in_weighted          <- read_and_clean("in_weighted.csv", gcp)
  cov_out_weighted         <- read_and_clean("out_weighted.csv", gcp)
  cov_positiveinfluence    <- read_and_clean("positiveinfluence.csv", gcp)
  cov_negativeinfluence    <- read_and_clean("negativeinfluence.csv", gcp)
  cov_positiveinf_weighted <- read_and_clean("positiveinf_weighted.csv", gcp)
  cov_negativeinf_weighted <- read_and_clean("negativeinf_weighted.csv", gcp)
  cov_distance             <- read_and_clean("distance_traveled_between_interactions.csv", gcp)
  cov_n_soc                <- read_and_clean("number_of_flies_in_soc_space.csv", gcp)
  cov_uniq_int             <- read_and_clean("unique_partners_met_interaction_space.csv", gcp)
  cov_uniq_soc             <- read_and_clean("unique_partners_met_social_space.csv", gcp)

  nodesAttr           <- make_nodes(flies)
  interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
  interaction_network <- link_events(x = interaction_network, change_event = interaction_data, nodes = nodesAttr)
  dependent           <- make_dependent_events(events = interaction_data, nodes = nodesAttr,
                                               default_network = interaction_network)

  nodesAttr <- link_events(nodesAttr, cov_in_degree, "popularity")
  nodesAttr <- link_events(nodesAttr, cov_out_degree, "activity")
  nodesAttr <- link_events(nodesAttr, cov_in_weighted, "popularity_weighted")
  nodesAttr <- link_events(nodesAttr, cov_out_weighted, "activity_weighted")
  nodesAttr <- link_events(nodesAttr, cov_positiveinfluence, "positive_influence")
  nodesAttr <- link_events(nodesAttr, cov_negativeinfluence, "negative_influence")
  nodesAttr <- link_events(nodesAttr, cov_positiveinf_weighted, "positive_inf_weighted")
  nodesAttr <- link_events(nodesAttr, cov_negativeinf_weighted, "negative_inf_weighted")
  nodesAttr <- link_events(nodesAttr, cov_distance, "distance_traveled_between_interactions")
  nodesAttr <- link_events(nodesAttr, cov_n_soc, "number_of_flies_in_soc_space")
  nodesAttr <- link_events(nodesAttr, cov_uniq_int, "unique_partners_met_interaction_space")
  nodesAttr <- link_events(nodesAttr, cov_uniq_soc, "unique_partners_met_social_space")

  dataDynam <- make_data(nodesAttr, interaction_network, dependent, flies)

  data_rate <- make_data_re(
    random_effects = NULL,
    fixed_effects  = rate_model_formula,
    model = "DyNAM", sub_model = "rate",
    data = dataDynam,
    control_preprocessing = set_preprocessing_opt(start_time = 0)
  )

  cat(sprintf("%s [PID %d] END   %2d/%d  %s  (%.1fs)\n", format(Sys.time(), "%F %T"),
              Sys.getpid(), i, nrow(files_to_process), group_name,
              as.numeric(difftime(Sys.time(), .t0, units = "secs"))), file = group_log, append = TRUE)

  list(group_name = group_name, treatment_name = treatment_name, group_index = i, data_object = data_rate)
}
stopCluster(cl); registerDoSEQ()

all_groups_data_list <- list(); all_groups_info_list <- list()
for (result in parallel_results) {
  if (inherits(result, "list") && !is.null(result$data_object)) {
    all_groups_data_list[[result$group_name]] <- result$data_object
    all_groups_info_list[[result$group_name]] <- data.frame(
      group = result$group_name, treatment = result$treatment_name, ixGroup = result$group_index)
  } else log_msg("GROUP FAILED: ", paste(utils::capture.output(print(result)), collapse = " "))
}
all_groups_info <- do.call(rbind, all_groups_info_list)

saveRDS(all_groups_data_list, file.path(out_path, "groups_data_rate.rds"))
saveRDS(all_groups_info,      file.path(out_path, "groups_info_rate.rds"))
log_msg("=== rate prep DONE ===  ", length(all_groups_data_list), "/", nrow(files_to_process),
        " groups saved to ", out_path)
log_msg("treatment split: ", paste(names(table(all_groups_info$treatment)),
        table(all_groups_info$treatment), sep = "=", collapse = "  "))
message("\nNext: Rscript rate_2_fit.r")
