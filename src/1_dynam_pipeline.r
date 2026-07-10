#!/usr/bin/env Rscript
# =============================================================================
# CHOICE DyNAM-RE  — fresh run on goldfish 1.8.6 / goldfish.latent 0.1.3
# Run it in ONE shot:   Rscript run_choice.r
# (do NOT paste line-by-line into the console)
#
# Ignore every old saved .rds — those are from the broken versions.
# This script does 3 things: (1) preprocess groups, (2) gather+scale,
# (3) compile Stan + sample, then saves everything to res/choice/.
# =============================================================================

# Develop versions of the packages used in this script:
# remotes::install_github("stocnet/goldfish@refactor/rate_prep")
# remotes::install_github('auzaheta/goldfish.latent@develop')

suppressPackageStartupMessages({
  library(goldfish.latent)
  library(goldfish)
  library(here)
  library(glue)
  library(cmdstanr)
  library(foreach)
  library(doParallel)
  library(parallel)
})
# ===== SET THESE ONCE =====
PROJECT <- "/home/milky/droso-pipe/6_dynam_models" # folder containing data/  (confirmed correct)
SRC_DIR <- file.path(PROJECT, "src") # folder containing utils.R + the .stan  <-- EDIT if `find` shows elsewhere

base_path <- PROJECT
edgelist_path <- file.path(base_path, "data", "edgelists")
covariances_path <- file.path(base_path, "data", "covariances")
results_path <- file.path(base_path, "res")
out_path <- file.path(results_path, "choice")
utils_file <- file.path(SRC_DIR, "utils.R")
stan_file <- system.file(
  "stan",
  "DNRE_Q1_choice_int_centered.stan",
  package = "goldfish.latent"
)
treatments_list <- c("CS_10D", "Cs_5DIZ", "CsCh")

if (!dir.exists(edgelist_path)) {
  stop(
    "No data/edgelists under PROJECT -> fix the PROJECT path at the top.\n  looked in: ",
    edgelist_path
  )
}
if (!file.exists(utils_file)) {
  stop("utils.R not found -> fix PROJECT.\n  looked in: ", utils_file)
}
source(utils_file) # read_and_clean() and other helpers
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# ---- fixed-effects formula -------------------------------------------------
choice_model_formula <- dependent ~
  indeg(interaction_network, weighted = FALSE, window = 96) +
  indeg(interaction_network, weighted = FALSE, window = 288) +
  indeg(interaction_network, weighted = FALSE, window = 864) +
  indeg(interaction_network, weighted = FALSE, window = 2592) +
  outdeg(interaction_network, weighted = FALSE, window = 96) +
  outdeg(interaction_network, weighted = FALSE, window = 288) +
  outdeg(interaction_network, weighted = FALSE, window = 864) +
  outdeg(interaction_network, weighted = FALSE, window = 2592) +
  trans(interaction_network, window = 96) +
  trans(interaction_network, window = 288) +
  trans(interaction_network, window = 864) +
  trans(interaction_network, window = 2592) +
  inertia(interaction_network, weighted = TRUE, window = 96) +
  inertia(interaction_network, weighted = TRUE, window = 288) +
  inertia(interaction_network, weighted = TRUE, window = 864) +
  inertia(interaction_network, weighted = TRUE, window = 2592) +
  inertia(interaction_network, weighted = FALSE, window = 96) +
  # window 288 (weighted=FALSE) is the RANDOM effect -> kept out of fixed effects
  inertia(interaction_network, weighted = FALSE, window = 864) +
  inertia(interaction_network, weighted = FALSE, window = 2592) +
  recip(interaction_network, weighted = FALSE, window = 96) +
  recip(interaction_network, weighted = FALSE, window = 288) +
  recip(interaction_network, weighted = FALSE, window = 864) +
  recip(interaction_network, weighted = FALSE, window = 2592) +
  alter(nodesAttr$number_of_flies_in_soc_space) +
  sim(nodesAttr$activity) +
  alter(flies$young) +
  alter(flies$old) +
  alter(flies$isolated)

# ===========================================================================
# STEP 1 — preprocess every group (parallel)
# ===========================================================================
all_files_list <- lapply(treatments_list, function(treatment) {
  folder_path <- file.path(edgelist_path, treatment)
  csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  if (length(csv_files) > 0) {
    data.frame(
      treatment_name = treatment,
      file_path = csv_files,
      group_name = tools::file_path_sans_ext(basename(csv_files))
    )
  } else {
    NULL
  }
})
files_to_process <- do.call(rbind, all_files_list)

n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

parallel_results <- foreach(
  i = seq_len(nrow(files_to_process)),
  .packages = c("goldfish", "goldfish.latent", "glue", "tools"),
  .errorhandling = "pass"
) %dopar%
  {
    treatment_name <- files_to_process$treatment_name[i]
    group_name <- files_to_process$group_name[i]
    file_path <- files_to_process$file_path[i]

    interaction_data <- read.csv(file_path)
    interaction_data <- interaction_data[, c(
      "time",
      "sender",
      "receiver",
      "increment"
    )]

    # getActors() was REMOVED in goldfish 1.8.x -> build the node frame directly.
    actors <- sort(unique(c(
      interaction_data$sender,
      interaction_data$receiver
    )))
    flies <- data.frame(
      label = as.character(actors),
      present = TRUE,
      stringsAsFactors = FALSE
    )
    for (col in c(
      "popularity",
      "activity",
      "popularity_weighted",
      "activity_weighted",
      "positive_influence",
      "negative_influence",
      "positive_inf_weighted",
      "negative_inf_weighted",
      "distance_traveled_between_interactions",
      "number_of_flies_in_soc_space",
      "unique_partners_met_interaction_space",
      "unique_partners_met_social_space",
      "young",
      "old",
      "isolated"
    )) {
      flies[[col]] <- 0
    }
    if (treatment_name == "CsCh") {
      flies$young <- 1
    }
    if (treatment_name == "CS_10D") {
      flies$old <- 1
    }
    if (treatment_name == "Cs_5DIZ") {
      flies$isolated <- 1
    }

    group_cov_path <- file.path(covariances_path, treatment_name, group_name)
    cov_in_degree <- read_and_clean("in_degree.csv", group_cov_path)
    cov_out_degree <- read_and_clean("out_degree.csv", group_cov_path)
    cov_in_weighted <- read_and_clean("in_weighted.csv", group_cov_path)
    cov_out_weighted <- read_and_clean("out_weighted.csv", group_cov_path)
    cov_positiveinfluence <- read_and_clean(
      "positiveinfluence.csv",
      group_cov_path
    )
    cov_negativeinfluence <- read_and_clean(
      "negativeinfluence.csv",
      group_cov_path
    )
    cov_positiveinf_weighted <- read_and_clean(
      "positiveinf_weighted.csv",
      group_cov_path
    )
    cov_negativeinf_weighted <- read_and_clean(
      "negativeinf_weighted.csv",
      group_cov_path
    )
    cov_distance <- read_and_clean(
      "distance_traveled_between_interactions.csv",
      group_cov_path
    )
    cov_n_soc <- read_and_clean(
      "number_of_flies_in_soc_space.csv",
      group_cov_path
    )
    cov_uniq_int <- read_and_clean(
      "unique_partners_met_interaction_space.csv",
      group_cov_path
    )
    cov_uniq_soc <- read_and_clean(
      "unique_partners_met_social_space.csv",
      group_cov_path
    )

    nodesAttr <- make_nodes(flies)
    interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
    interaction_network <- link_events(
      x = interaction_network,
      change_event = interaction_data,
      nodes = nodesAttr
    )
    dependent <- make_dependent_events(
      events = interaction_data,
      nodes = nodesAttr,
      default_network = interaction_network
    )

    nodesAttr <- link_events(nodesAttr, cov_in_degree, "popularity")
    nodesAttr <- link_events(nodesAttr, cov_out_degree, "activity")
    nodesAttr <- link_events(nodesAttr, cov_in_weighted, "popularity_weighted")
    nodesAttr <- link_events(nodesAttr, cov_out_weighted, "activity_weighted")
    nodesAttr <- link_events(
      nodesAttr,
      cov_positiveinfluence,
      "positive_influence"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_negativeinfluence,
      "negative_influence"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_positiveinf_weighted,
      "positive_inf_weighted"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_negativeinf_weighted,
      "negative_inf_weighted"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_distance,
      "distance_traveled_between_interactions"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_n_soc,
      "number_of_flies_in_soc_space"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_uniq_int,
      "unique_partners_met_interaction_space"
    )
    nodesAttr <- link_events(
      nodesAttr,
      cov_uniq_soc,
      "unique_partners_met_social_space"
    )

    dataDynam <- make_data(nodesAttr, interaction_network, dependent, flies)

    data_choice <- make_data_re(
      random_effects = list(
        inertia(interaction_network, weighted = FALSE, window = 288) ~ 1
      ),
      fixed_effects = choice_model_formula,
      model = "DyNAM",
      sub_model = "choice",
      data = dataDynam
    )

    list(
      group_name = group_name,
      treatment_name = treatment_name,
      group_index = i,
      data_object = data_choice
    )
  }
stopCluster(cl)
registerDoSEQ()

all_groups_data_list <- list()
all_groups_info_list <- list()
for (result in parallel_results) {
  if (inherits(result, "list") && !is.null(result$data_object)) {
    all_groups_data_list[[result$group_name]] <- result$data_object
    all_groups_info_list[[result$group_name]] <- data.frame(
      group = result$group_name,
      treatment = result$treatment_name,
      ixGroup = result$group_index
    )
  } else {
    message(
      "Group failed: ",
      paste(utils::capture.output(print(result)), collapse = " ")
    )
  }
}
all_groups_info <- do.call(rbind, all_groups_info_list)

saveRDS(all_groups_data_list, file.path(out_path, "groups_data.rds"))
saveRDS(all_groups_info, file.path(out_path, "groups_info.rds"))
message("STEP 1 done: ", length(all_groups_data_list), " groups preprocessed.")

# ===========================================================================
# STEP 2 — gather + build interaction design matrix + scale
# ===========================================================================
data_stan_re <- gather_groups(
  groups_info = all_groups_info,
  data = all_groups_data_list
)

# guard: fail loud if effect names shifted (would silently break the slicing below)
.cn <- colnames(data_stan_re$data_stan$X_choice)
cat_interaction_names <- grep("alter_fliesOf", .cn, value = TRUE, fixed = TRUE)
cat_interaction_pos <- grep("alter_fliesOf", .cn, value = FALSE, fixed = TRUE)
stopifnot(
  "sim_nodesAttrOfpopularity" %in% .cn,
  any(cat_interaction_names %in% .cn)
)

# transform the distance covariate (log1p) due to heavy tail
var_distance <- grep("distance", .cn, value = TRUE)
data_stan_re$data_stan$X_choice[, var_distance] <-
  log1p(data_stan_re$data_stan$X_choice[, var_distance])

interaction_all <- data_stan_re$data_stan$X_choice[, cat_interaction_names] %*%
  seq_along(cat_interaction_names)

data_stan_re$data_stan$X_choice <- data_stan_re$data_stan$X_choice[,
  -cat_interaction_pos,
  drop = FALSE
]

names(data_stan_re$data_stan)[
  names(data_stan_re$data_stan) == "X_choice"
] <- "X_choice_raw"

data_stan_re$data_stan$grain_size <- 4
data_stan_re$data_stan$V1_choice <-
  grep(colnames(data_stan_re$data_stan$Z_choice), .cn, fixed = TRUE)
data_stan_re$data_stan$P_choice <- ncol(data_stan_re$data_stan$X_choice_raw)

# interaction data
data_stan_re$data_stan$send_int <- rep(1:3, each = 20)
sg <- data_stan_re$data_stan$start_group
data_stan_re$data_stan$interaction <- c(
  rep(1, sg[20]),
  rep(2, sg[40] - sg[20]),
  rep(3, data_stan_re$data_stan$T_choice - sg[40])
)
data_stan_re$data_stan$C <- 3

# clean up
data_stan_re$data_stan$Z_choice <- NULL

dim_x_choice <- dim(data_stan_re$data_stan$X_choice_raw)

message(
  "STEP 2 done: design matrix ",
  dim_x_choice[1],
  " x ",
  dim_x_choice[2]
)

# ===========================================================================
# STEP 3 — compile Stan + sample
# ===========================================================================
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = data_stan_re$data_stan,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 3,
  iter_warmup = 100,
  iter_sampling = 100
)

# ---- save everything needed for later analysis (no reload guesswork) ------
fit$save_object(file.path(out_path, "fit_choice.rds")) # self-contained
saveRDS(fit$draws(), file.path(out_path, "draws_choice.rds")) # posterior array
saveRDS(data_stan_re, file.path(out_path, "data_stan_re_choice.rds")) # <- SAVED THIS TIME
message("DONE. Choice fit + draws + data_stan_re saved to: ", out_path)
