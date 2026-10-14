#!/usr/bin/env Rscript
# =============================================================================
# Validate the statistic engine against goldfish's OWN output.
# Replays the observed events of one group and compares the computed choice
# statistics to X_choice_raw column-for-column. If this does not match, the
# simulator built on top of it would be silently wrong.
#   Rscript src/5b_validate_engine.r [group_index]
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
source("src/utils.R")
source("src/5a_stat_engine.r")

a <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(a)) a <- 1L

gi  <- readRDS("res/choice/groups_info.rds")
ds  <- readRDS("res/choice/data_stan_re_choice.rds")$data_stan
grp <- as.character(gi$group[a]); trt <- as.character(gi$treatment[a])
cat("validating group", a, grp, "(", trt, ")\n")

ev <- read.csv(file.path("data/edgelists", trt, paste0(grp, ".csv")))
ev <- ev[order(ev$time), c("time", "sender", "receiver", "increment")]
labels <- as.character(sort(unique(c(ev$sender, ev$receiver))))
n <- length(labels)
si <- match(as.character(ev$sender), labels)
ri <- match(as.character(ev$receiver), labels)
E  <- nrow(ev)
cat("events:", E, " actors:", n, "\n")

# rows of X for this group
ev_lo <- ds$start_group[a]
ev_hi <- if (a < length(ds$start_group)) ds$start_group[a + 1] - 1L else ds$T_choice
stopifnot(ev_hi - ev_lo + 1L == E)
row_lo <- ds$start_choice[ev_lo]; row_hi <- ds$end_choice[ev_hi]
Xgf <- ds$X_choice_raw[row_lo:row_hi, , drop = FALSE]
cat("goldfish rows for this group:", nrow(Xgf), " (expect", E * (n - 1), ")\n")

# covariates
gcp <- file.path("data/covariances", trt, grp)
cov_files <- c(activity = "out_degree.csv",
               number_of_flies_in_soc_space = "number_of_flies_in_soc_space.csv")
cs <- lapply(cov_files, function(f) new_cov_state(read_and_clean(f, gcp), labels))

spec <- choice_spec()
WS <- new_window_state(CHOICE_WINDOWS, n)
events <- list(time = ev$time, s = si, r = ri, inc = ev$increment, n = E)

Xme <- matrix(0, E * (n - 1), nrow(spec))
pos <- 1L
for (e in seq_len(E)) {
  t <- ev$time[e]
  for (k in seq_along(WS)) WS[[k]] <- expire(WS[[k]], events, t)
  for (nm in names(cs))    cs[[nm]] <- advance_cov(cs[[nm]], t)
  cov_cur <- lapply(cs, `[[`, "cur")

  i   <- si[e]
  alt <- setdiff(seq_len(n), i)
  Xme[pos:(pos + n - 2), ] <- choice_stats(i, alt, WS, cov_cur, spec)
  pos <- pos + n - 1L

  for (k in seq_along(WS)) WS[[k]] <- add_event(WS[[k]], si[e], ri[e], ev$increment[e])
}

cn <- colnames(ds$X_choice_raw)
cat("\n", strrep("=", 78), "\nCOLUMN-BY-COLUMN COMPARISON vs goldfish\n",
    strrep("=", 78), "\n", sep = "")
res <- data.frame(col = seq_along(cn), name = cn,
                  max_abs_diff = NA_real_, exact = NA, stringsAsFactors = FALSE)
for (p in seq_along(cn)) {
  d <- max(abs(Xme[, p] - Xgf[, p]))
  res$max_abs_diff[p] <- d
  res$exact[p] <- d < 1e-8
}
print(res, row.names = FALSE, digits = 4)
cat("\nMATCHING:", sum(res$exact), "/", nrow(res), "columns\n")
if (all(res$exact)) {
  cat("\nENGINE VALIDATED — statistics reproduce goldfish exactly.\n")
} else {
  cat("\nMISMATCH in:", paste(res$name[!res$exact], collapse = ", "), "\n")
}
saveRDS(list(Xme = Xme, Xgf = Xgf, res = res), "/tmp/claude-1030/-home-milky-droso-pipe-6-dynam-models/80b702ca-deda-4e81-b784-ffcc2ba238a6/scratchpad/validate.rds")
