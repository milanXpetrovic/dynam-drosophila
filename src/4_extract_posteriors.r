#!/usr/bin/env Rscript
# =============================================================================
# Extract posterior distributions for BOTH models to plain CSV.
#   Rscript src/4_extract_posteriors.r
# Writes to res/posterior/:
#   <model>_summary.csv   one row per parameter: median/mean/sd/90% CrI/Rhat/ESS
#   <model>_draws.csv.gz  every draw of every parameter (wide: 1 row = 1 draw)
# Parameters are labelled with predictor / treatment / group names, not indices.
# =============================================================================

suppressPackageStartupMessages({
  library(posterior)
})

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
res     <- file.path(PROJECT, "res")
outdir  <- file.path(res, "posterior")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

TREATMENTS <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")

# 90% CrI (q5/q95) and 95% CrI (q2.5/q97.5)
q2 <- function(...) posterior::quantile2(..., probs = c(0.025, 0.05, 0.95, 0.975))

# index i and j out of "name[i,j]" / "name[i]" ------------------------------
idx <- function(v, which_one) {
  m <- regmatches(v, regexec("\\[([0-9]+)(?:,([0-9]+))?\\]", v))
  as.integer(vapply(m, function(z) if (length(z) >= which_one + 1) z[which_one + 1] else NA_character_, ""))
}

extract_model <- function(model) {
  fit_file <- file.path(res, model, sprintf("fit_%s.rds", model))
  dsr_file <- file.path(res, model, sprintf("data_stan_re_%s.rds", model))
  gi_file  <- file.path(res, model,
                        if (model == "choice") "groups_info.rds" else "groups_info_rate.rds")
  for (f in c(fit_file, dsr_file, gi_file))
    if (!file.exists(f)) stop("missing input: ", f)

  fit <- readRDS(fit_file)
  dsr <- readRDS(dsr_file)
  gi  <- readRDS(gi_file)

  # predictor names, minus the random-effect column where there is one
  if (model == "choice") {
    xcols  <- colnames(dsr$data_stan$X_choice_raw)
    re_col <- xcols[dsr$data_stan$V1_choice]
    pnames <- xcols[-dsr$data_stan$V1_choice]
    beta   <- "beta_choice"
  } else {
    pnames <- colnames(dsr$data_stan$X_rate_raw)
    re_col <- "(intercept)"
    beta   <- "beta_rate"
  }
  groups <- as.character(gi$group)

  s <- fit$summary(NULL, mean, median, sd, mad, q2, default_convergence_measures())
  fam <- sub("\\[.*", "", s$variable)
  i1  <- idx(s$variable, 1)
  i2  <- idx(s$variable, 2)

  # human-readable labels ---------------------------------------------------
  s$family    <- fam
  s$treatment <- NA_character_
  s$group     <- NA_character_
  s$predictor <- NA_character_
  s$meaning   <- NA_character_

  is_beta <- fam == beta
  s$treatment[is_beta] <- TREATMENTS[i1[is_beta]]
  s$predictor[is_beta] <- pnames[i2[is_beta]]
  s$meaning[is_beta]   <- "fixed effect (per treatment)"

  is_gamma <- fam == "gamma"
  s$group[is_gamma]   <- groups[i1[is_gamma]]
  s$meaning[is_gamma] <- if (model == "choice")
    sprintf("group random effect on %s", re_col) else "group random intercept"

  is_mu <- fam == "mu"
  s$treatment[is_mu] <- TREATMENTS[i1[is_mu]]
  s$meaning[is_mu]   <- sprintf("mean of group random effect on %s", re_col)

  is_alpha <- fam == "alpha"
  s$treatment[is_alpha] <- TREATMENTS[i1[is_alpha]]
  s$meaning[is_alpha]   <- "intercept (per treatment)"

  s$meaning[fam == "sigma"] <- "SD of the group random effects"
  s$meaning[fam == "lp__"]  <- "log posterior density"
  s$meaning[grepl("_std$", fam)] <- "non-centred (standardised) parameterisation"

  front <- c("variable", "family", "meaning", "treatment", "predictor", "group")
  s <- s[, c(front, setdiff(names(s), front))]

  sum_file <- file.path(outdir, sprintf("%s_summary.csv", model))
  write.csv(s, sum_file, row.names = FALSE)

  # full draws --------------------------------------------------------------
  dd <- as_draws_df(fit$draws())
  draw_file <- file.path(outdir, sprintf("%s_draws.csv.gz", model))
  con <- gzfile(draw_file, "w"); write.csv(as.data.frame(dd), con, row.names = FALSE); close(con)

  message(sprintf("%-6s -> %s (%d params) + %s (%d draws x %d cols)",
                  model, basename(sum_file), nrow(s), basename(draw_file),
                  nrow(dd), ncol(dd)))
  invisible(s)
}

for (m in c("choice", "rate")) extract_model(m)

message("\nWrote to ", outdir)
print(data.frame(file = list.files(outdir),
                 size = sprintf("%.1f MB", file.size(file.path(outdir, list.files(outdir))) / 1e6)))
