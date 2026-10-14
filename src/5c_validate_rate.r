#!/usr/bin/env Rscript
# =============================================================================
# Validate the RATE statistics against goldfish's X_rate_raw.
# The rate design discretises time into unit steps: interval k covers [k-1, k),
# with one row per actor. Statistics are the state at the start of the interval.
#   Rscript src/5c_validate_rate.r [group_index]
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
source("src/utils.R")
source("src/5a_stat_engine.r")

a <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(a)) a <- 1L

gi  <- readRDS("res/rate/groups_info_rate.rds")
ds  <- readRDS("res/rate/data_stan_re_rate.rds")$data_stan
grp <- as.character(gi$group[a]); trt <- as.character(gi$treatment[a])
cat("validating RATE for group", a, grp, "(", trt, ")\n")

ev <- read.csv(file.path("data/edgelists", trt, paste0(grp, ".csv")))
ev <- ev[order(ev$time), c("time", "sender", "receiver", "increment")]
labels <- as.character(sort(unique(c(ev$sender, ev$receiver))))
n <- length(labels)
si <- match(as.character(ev$sender), labels)
ri <- match(as.character(ev$receiver), labels)
E  <- nrow(ev)

iv_lo <- ds$start_group[a]
iv_hi <- if (a < length(ds$start_group)) ds$start_group[a + 1] - 1L else ds$T_rate
n_iv  <- iv_hi - iv_lo + 1L
row_lo <- (iv_lo - 1L) * n + 1L
row_hi <- iv_hi * n
Xgf <- ds$X_rate_raw[row_lo:row_hi, , drop = FALSE]
cat("intervals:", n_iv, " actors:", n, " rows:", nrow(Xgf), "\n")

# covariates used by the rate model (log1p on distance, as in 3_rate_fit.r)
cov_files <- c(
  positive_influence    = "positiveinfluence.csv",
  negative_influence    = "negativeinfluence.csv",
  positive_inf_weighted = "positiveinf_weighted.csv",
  negative_inf_weighted = "negativeinf_weighted.csv",
  distance_traveled_between_interactions = "distance_traveled_between_interactions.csv",
  number_of_flies_in_soc_space = "number_of_flies_in_soc_space.csv",
  unique_partners_met_social_space = "unique_partners_met_social_space.csv"
)
gcp <- file.path("data/covariances", trt, grp)
cs <- lapply(names(cov_files), function(nm) {
  df <- read_and_clean(cov_files[[nm]], gcp)
  if (nm == "distance_traveled_between_interactions") df$replace <- log1p(df$replace)
  new_cov_state(df, labels)
})
names(cs) <- names(cov_files)

spec <- rate_spec()
WS <- new_window_state(CHOICE_WINDOWS, n)
events <- list(time = ev$time, s = si, r = ri, inc = ev$increment, n = E)

Xme <- matrix(0, n_iv * n, nrow(spec))
ei <- 1L
for (k in seq_len(n_iv)) {
  t <- k - 1L                                  # interval k covers [k-1, k)
  for (w in seq_along(WS)) WS[[w]] <- expire(WS[[w]], events, t)
  for (nm in names(cs))    cs[[nm]] <- advance_cov(cs[[nm]], t)
  cov_cur <- lapply(cs, `[[`, "cur")

  Xme[((k - 1L) * n + 1L):(k * n), ] <- rate_stats(WS, cov_cur, spec, n)

  while (ei <= E && ev$time[ei] == t) {         # apply events occurring at t
    for (w in seq_along(WS)) WS[[w]] <- add_event(WS[[w]], si[ei], ri[ei], ev$increment[ei])
    ei <- ei + 1L
  }
}

cn <- colnames(ds$X_rate_raw)
cat("\n", strrep("=", 70), "\nRATE COLUMN COMPARISON vs goldfish\n", strrep("=", 70), "\n", sep = "")
res <- data.frame(col = seq_along(cn), name = cn, max_abs_diff = NA_real_,
                  exact = NA, stringsAsFactors = FALSE)
for (p in seq_along(cn)) {
  d <- max(abs(Xme[, p] - Xgf[, p]))
  res$max_abs_diff[p] <- d
  res$exact[p] <- d < 1e-8
}
print(res, row.names = FALSE, digits = 4)
cat("\nMATCHING:", sum(res$exact), "/", nrow(res), "columns\n")
if (all(res$exact)) cat("\nRATE ENGINE VALIDATED\n") else
  cat("\nMISMATCH in:", paste(res$name[!res$exact], collapse = ", "), "\n")
