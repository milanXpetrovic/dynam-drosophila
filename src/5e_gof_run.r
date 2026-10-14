#!/usr/bin/env Rscript
# =============================================================================
# Posterior predictive check by forward simulation.
#   nohup Rscript src/5e_gof_run.r [n_draws] > logs/5_gof.log 2>&1 &
#
# For every group (60) and every posterior draw (default 100) simulate one
# event sequence over the observed time horizon, then compute GoF statistics.
# Results are aggregated per treatment across the groups within it.
#
# Built on the VALIDATED statistic engine (5a): choice 26/26 and rate 9/9
# columns reproduce goldfish exactly.
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
suppressPackageStartupMessages({
  library(posterior); library(foreach); library(doParallel); library(parallel)
})
source("src/utils.R"); source("src/5a_stat_engine.r"); source("src/5d_simulate.r")

N_DRAWS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(N_DRAWS)) N_DRAWS <- 100L
OUT <- file.path("res", "gof"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TREATMENTS <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")

log_msg <- function(...) message(format(Sys.time(), "%F %T"), "  ", ...)
log_msg("=== GoF START | draws = ", N_DRAWS, " ===")

gi     <- readRDS("res/choice/groups_info.rds")
cc_all <- readRDS("res/choice/data_stan_re_choice.rds")$data_stan$send_int
dc <- as_draws_matrix(readRDS("res/choice/fit_choice.rds")$draws())
dr <- as_draws_matrix(readRDS("res/rate/fit_rate.rds")$draws())
# thin each posterior evenly to N_DRAWS (the two models were fitted separately,
# so draw s of one carries no dependence on draw s of the other)
ix_c <- round(seq(1, nrow(dc), length.out = N_DRAWS))
ix_r <- round(seq(1, nrow(dr), length.out = N_DRAWS))
DR <- list(choice = dc[ix_c, , drop = FALSE], rate = dr[ix_r, , drop = FALSE])
rm(dc, dr); gc()

# ---- GoF statistics ---------------------------------------------------------
gini <- function(x) {
  x <- sort(x); n <- length(x); s <- sum(x)
  if (s == 0) return(0)
  sum((2 * seq_len(n) - n - 1) * x) / (n * s)
}

gof_stats <- function(ev, labels) {
  N <- length(labels)
  base <- network_stats(ev, labels)
  if (nrow(ev) < 3) {
    return(c(base[c("recip_ratio", "w_density", "n_events")],
             list(indeg_sd = 0, outdeg_sd = 0, indeg_gini = 0, outdeg_gini = 0,
                  tie_w_mean = 0, tie_w_sd = 0, tie_w_max = 0,
                  iet_median = NA, iet_cv = NA, burstiness = NA,
                  repeat_frac = NA, two_paths = 0, transitivity = NA)))
  }
  # temporal
  iet <- diff(sort(ev$time))
  burst <- (sd(iet) - mean(iet)) / (sd(iet) + mean(iet))
  # repetition: consecutive events on the same ordered dyad
  dy <- paste(ev$sender, ev$receiver, sep = "->")
  rep_frac <- mean(dy[-1] == dy[-length(dy)])
  # aggregated binary network -> two-paths and global transitivity
  B <- matrix(0, N, N)
  B[cbind(match(ev$sender, labels), match(ev$receiver, labels))] <- 1
  diag(B) <- 0
  TP <- B %*% B
  two_paths <- sum(TP) - sum(diag(TP))
  closed <- sum(TP * B)
  list(recip_ratio = base$recip_ratio, w_density = base$w_density,
       n_events = base$n_events,
       indeg_sd = sd(base$indeg), outdeg_sd = sd(base$outdeg),
       indeg_gini = gini(base$indeg), outdeg_gini = gini(base$outdeg),
       tie_w_mean = mean(base$tie_weights), tie_w_sd = sd(base$tie_weights),
       tie_w_max = max(base$tie_weights),
       iet_median = median(iet), iet_cv = sd(iet) / mean(iet), burstiness = burst,
       repeat_frac = rep_frac, two_paths = two_paths,
       transitivity = if (two_paths > 0) closed / two_paths else NA)
}

# distribution summaries, pooled later per treatment
QS <- seq(0, 1, by = 0.1)
dist_summary <- function(ev, labels) {
  b <- network_stats(ev, labels)
  list(indeg = b$indeg, outdeg = b$outdeg, tie_weights = b$tie_weights,
       iet = if (nrow(ev) > 1) diff(sort(ev$time)) else numeric(0))
}

# ---- run --------------------------------------------------------------------
n_cores <- max(1L, parallel::detectCores() - 2L)
cl <- makeCluster(n_cores); registerDoParallel(cl)
log_msg("cores: ", n_cores, " | groups: ", nrow(gi), " | sims: ", nrow(gi) * N_DRAWS)

results <- foreach(
  a = seq_len(nrow(gi)),
  .packages = c("posterior"),
  .export = c("build_group", "simulate_group", "network_stats", "gof_stats",
              "dist_summary", "make_params", "gini", "choice_spec", "rate_spec",
              "new_window_state", "expire", "add_event", "new_cov_state",
              "advance_cov", "choice_stats", "rate_stats", "bin",
              "CHOICE_WINDOWS", "read_and_clean", "DR", "gi", "cc_all", "N_DRAWS"),
  .errorhandling = "pass"
) %dopar% {
  G  <- build_group(a, gi)
  cc <- cc_all[a]
  sc <- choice_spec(); sr <- rate_spec()

  obs_stats <- gof_stats(G$obs, G$labels)
  obs_dist  <- dist_summary(G$obs, G$labels)

  sim_stats <- vector("list", N_DRAWS)
  sim_dist  <- vector("list", N_DRAWS)
  for (s in seq_len(N_DRAWS)) {
    set.seed(1000L * a + s)                       # reproducible per (group, draw)
    p <- make_params(a, cc, s, DR)
    sim <- simulate_group(G, p$rate, p$choice, spec_rate = sr, spec_choice = sc)
    sim_stats[[s]] <- gof_stats(sim, G$labels)
    sim_dist[[s]]  <- dist_summary(sim, G$labels)
  }
  list(index = a, group = G$group, treatment = G$treatment, cc = cc,
       n_obs = G$n_obs, obs = obs_stats, obs_dist = obs_dist,
       sim = sim_stats, sim_dist = sim_dist)
}
stopCluster(cl); registerDoSEQ()

ok <- vapply(results, function(r) is.list(r) && !is.null(r$obs), TRUE)
if (any(!ok)) {
  log_msg("FAILED groups: ", sum(!ok))
  for (r in results[!ok]) log_msg("  ", paste(utils::capture.output(print(r)), collapse = " "))
}
results <- results[ok]
log_msg("simulated ", length(results), "/", nrow(gi), " groups")

saveRDS(results, file.path(OUT, "gof_raw.rds"))

# ---- tidy scalar table ------------------------------------------------------
stat_names <- names(results[[1]]$obs)
rows <- do.call(rbind, lapply(results, function(r) {
  sim <- do.call(rbind, lapply(seq_along(r$sim), function(s)
    data.frame(draw = s, statistic = stat_names,
               value = unlist(r$sim[[s]])[stat_names], row.names = NULL)))
  sim$group <- r$group; sim$treatment <- TREATMENTS[r$cc]; sim$kind <- "simulated"
  obs <- data.frame(draw = NA_integer_, statistic = stat_names,
                    value = unlist(r$obs)[stat_names], row.names = NULL,
                    group = r$group, treatment = TREATMENTS[r$cc], kind = "observed")
  rbind(sim, obs)
}))
write.csv(rows, file.path(OUT, "gof_scalars.csv"), row.names = FALSE)

# ---- per-treatment aggregation ---------------------------------------------
# For each draw, average the statistic over the groups in a treatment -> a
# posterior predictive distribution for the treatment-level mean. Compare with
# the observed treatment-level mean and report a Bayesian p-value.
agg <- do.call(rbind, lapply(split(rows, rows$treatment), function(tr) {
  sim <- tr[tr$kind == "simulated", ]
  obs <- tr[tr$kind == "observed", ]
  do.call(rbind, lapply(unique(sim$statistic), function(st) {
    ss <- sim[sim$statistic == st, ]
    per_draw <- tapply(ss$value, ss$draw, mean, na.rm = TRUE)
    o <- mean(obs$value[obs$statistic == st], na.rm = TRUE)
    data.frame(treatment = tr$treatment[1], statistic = st,
               observed = o,
               sim_mean = mean(per_draw, na.rm = TRUE),
               sim_q05 = unname(quantile(per_draw, .05, na.rm = TRUE)),
               sim_q95 = unname(quantile(per_draw, .95, na.rm = TRUE)),
               ppp = mean(per_draw >= o, na.rm = TRUE),
               n_groups = length(unique(ss$group)))
  }))
}))
write.csv(agg, file.path(OUT, "gof_by_treatment.csv"), row.names = FALSE)

log_msg("=== GoF DONE === wrote gof_raw.rds, gof_scalars.csv, gof_by_treatment.csv to ", OUT)
print(agg[order(agg$treatment, agg$statistic), ], row.names = FALSE, digits = 3)
