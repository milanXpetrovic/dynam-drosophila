#!/usr/bin/env Rscript
# =============================================================================
# RATE stage 2 of 2 — FIT.  Run AFTER rate_1_prep.r (re-runnable):
#   cd /home/milky/droso-pipe/6_dynam_models
#   nohup Rscript rate_2_fit.r > res/rate/fit.log 2>&1 &
#   tail -f res/rate/fit.log
# Reloads groups_data_rate.rds from disk, builds the Stan data, samples.
# =============================================================================

suppressPackageStartupMessages({
  library(goldfish.latent); library(goldfish); library(glue); library(cmdstanr)
})

# ===== SET THESE ONCE (must match rate_1_prep.r) =====
PROJECT   <- "/home/milky/droso-pipe/6_dynam_models"
SRC_DIR   <- file.path(PROJECT, "src")
RATE_STAN <- file.path(SRC_DIR, "DNRE_Q0_rate_int_centered.stan")

base_path       <- PROJECT
out_path        <- file.path(base_path, "res", "rate")
utils_file      <- file.path(SRC_DIR, "utils.R")
treatments_list <- c("CS_10D", "Cs_5DIZ", "CsCh")

if (!file.exists(RATE_STAN)) stop("rate Stan not found -> fix RATE_STAN.\n  ", RATE_STAN)
if (file.exists(utils_file)) source(utils_file)

run_log <- file.path(out_path, format(Sys.time(), "fit_%Y%m%d_%H%M%S.log"))
log_msg <- function(...) { line <- paste0(format(Sys.time(), "%F %T"), "  ", ...)
  cat(line, "\n", sep = "", file = run_log, append = TRUE); message(line) }
log_msg("=== rate fit START ===")

# ---- reload STEP 1 output from disk ----------------------------------------
gd <- file.path(out_path, "groups_data_rate.rds")
gi <- file.path(out_path, "groups_info_rate.rds")
if (!file.exists(gd) || !file.exists(gi))
  stop("STEP 1 output not found under ", out_path, ".\n",
       "Run rate_1_prep.r first — it creates groups_data_rate.rds.")
all_groups_data_list <- readRDS(gd)
all_groups_info      <- readRDS(gi)
log_msg("reloaded ", length(all_groups_data_list), " groups from disk")

# ---- gather + build the custom-rate Stan data ------------------------------
data_stan_re <- gather_groups(groups_info = all_groups_info, data = all_groups_data_list)

old_X_rate <- data_stan_re$data_stan$X_rate
cn <- colnames(old_X_rate)
log_msg("X_rate columns (", length(cn), "):"); print(cn)

var_distance <- grep("distance", cn, value = TRUE)
old_X_rate[, var_distance] <- log1p(old_X_rate[, var_distance])

# robust removal: intercept, indeg/outdeg long windows (keep 96), node_trans,
# and the 3 treatment dummies (they become the categorical `interaction`).
pos_remove <-
  grepl("Intercept", cn) |
  grepl("^(indeg|outdeg)_interaction_network_(W_)?(288|864|2592)$", cn) |
  grepl("^node_?trans", cn) |
  cn %in% c("ego_fliesOfyoung", "ego_fliesOfold", "ego_fliesOfisolated")
stopifnot(
  any(grepl("Intercept", cn)),
  all(c("ego_fliesOfyoung", "ego_fliesOfold", "ego_fliesOfisolated") %in% cn),
  sum(grepl("^(indeg|outdeg)_interaction_network_(W_)?(288|864|2592)$", cn)) == 6
)
log_msg("removing ", sum(pos_remove), " cols, keeping ", sum(!pos_remove),
        " -> [", paste(cn[pos_remove], collapse = ", "), "]")

data_stan_re$data_stan$grain_size <- 1621   # NOTE: write_json(grain_size=) overrides in the JSON
data_stan_re$data_stan$Z_rate     <- NULL
data_stan_re$data_stan$X_rate_raw <- old_X_rate[, !pos_remove]
data_stan_re$data_stan$X_rate     <- NULL
data_stan_re$data_stan$P_rate     <- sum(!pos_remove)

# treatment-category structure (your design: 60 arenas, 20 per treatment, ordered)
log_msg("treatment split: ", paste(names(table(all_groups_info$treatment)),
        table(all_groups_info$treatment), sep = "=", collapse = "  "))
stopifnot(
  nrow(all_groups_info) == 60,
  identical(as.character(all_groups_info$treatment), rep(treatments_list, each = 20))
)
data_stan_re$data_stan$send_int <- rep(1:3, each = 20)
sg <- data_stan_re$data_stan$start_group
data_stan_re$data_stan$interaction <- c(
  rep(1, sg[20]),
  rep(2, sg[40] - sg[20]),
  rep(3, data_stan_re$data_stan$T_rate - sg[40])
)
data_stan_re$data_stan$C <- 3
data_stan_re$data_stan$sub_model <- NULL
attr(data_stan_re, "sub_model") <- "rate"

json_file <- goldfish.latent::write_json(data_stan = data_stan_re,
                        file_name = file.path(out_path, "rate.json"),
                        n_chunks = 100, grain_size = 10)
saveRDS(data_stan_re, file.path(out_path, "data_stan_re_rate.rds"))
log_msg("built Stan data: P_rate=", data_stan_re$data_stan$P_rate,
        "  T_rate=", data_stan_re$data_stan$T_rate)

# ---- compile + sample ------------------------------------------------------
log_msg("compiling ", RATE_STAN)
mod_rate <- cmdstan_model(RATE_STAN, cpp_options = list(stan_threads = TRUE))

csv_dir <- file.path(out_path, "stan_csv"); dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
log_msg("sampling started (the long one) — chain CSVs in ", csv_dir)

fit <- mod_rate$sample(
  data = json_file, chains = 4, parallel_chains = 4, threads_per_chain = 11,
  iter_warmup = 500, iter_sampling = 500, refresh = 50, output_dir = csv_dir
)

log_msg("sampler wall-times (s):")
cat(utils::capture.output(print(fit$time())), sep = "\n", file = run_log, append = TRUE)
fit$save_object(file.path(out_path, "fit_rate.rds"))
saveRDS(fit$draws(), file.path(out_path, "draws_rate.rds"))
log_msg("=== rate fit DONE ===  saved to ", out_path)