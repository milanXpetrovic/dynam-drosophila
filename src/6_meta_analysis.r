#!/usr/bin/env Rscript
# =============================================================================
# Random-effects meta-analysis of the per-group frequentist DyNAM fits, run
# separately within each treatment, to identify which effects vary most across
# groups -- i.e. the best CANDIDATES FOR GROUP-LEVEL RANDOM EFFECTS.
#
#   Rscript src/6_meta_analysis.r
#
# Input : old_res/old_results/{choice_model,rate_model}/<treatment>/*.json
#         (per-group estimates + standard errors from the frequentist fits)
# Output: res/meta/meta_<model>.csv        one row per treatment x effect
#         res/meta/meta_group_estimates.csv the per-group inputs
#         res/meta/plots/*.png              forest plots per effect
#
# tau^2 = between-group variance on the coefficient scale -> directly the prior
#         scale of a candidate random effect.
# I^2   = share of total variability that is between-group rather than sampling
#         error. High I^2 => a group-level random effect is warranted.
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
suppressPackageStartupMessages({ library(metafor); library(jsonlite) })

OUT <- file.path("res", "meta"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ROOT <- file.path("old_res", "old_results")
TREATMENT_LAB <- c(CS_10D = "CS_10D (old)", Cs_5DIZ = "Cs_5DIZ (isolated)", CsCh = "CsCh (young)")

# effect labels, in the order the parameters appear in each fit
EFFECTS <- list(
  choice_model = c("inertia", "recip_w288", "indeg_W_w864", "outdeg_W_w864",
                   "trans_w864", "sim_activity", "sim_popularity"),
  rate_model   = c("intercept", "indeg_W_w864", "outdeg_W_w864",
                   "ego_bursty", "ego_influence_pos", "ego_influence_neg")
)

DROP_NONCONVERGED <- TRUE

# ---- read every per-group fit ----------------------------------------------
read_fits <- function(model) {
  files <- list.files(file.path(ROOT, model), pattern = "\\.json$",
                      recursive = TRUE, full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    d <- jsonlite::fromJSON(f)
    est <- as.numeric(d$parameters); se <- as.numeric(d$standardErrors)
    nm <- EFFECTS[[model]]
    if (length(est) != length(nm))
      stop("unexpected parameter count in ", f, ": ", length(est), " vs ", length(nm))
    data.frame(
      model = model,
      treatment = TREATMENT_LAB[[basename(dirname(f))]],
      group = tools::file_path_sans_ext(basename(f)),
      effect = nm, estimate = est, se = se,
      converged = isTRUE(d$convergence$isConverged[1]),
      max_abs_score = as.numeric(d$convergence$maxAbsScore[1]),
      n_events = as.numeric(d$nEvents[1]),
      stringsAsFactors = FALSE
    )
  }))
}

fits <- rbind(read_fits("choice_model"), read_fits("rate_model"))
message("read ", nrow(fits), " group x effect estimates from ",
        length(unique(fits$group)) , " group files")

nc <- unique(fits[!fits$converged, c("model", "treatment", "group", "max_abs_score")])
if (nrow(nc)) {
  message("\nNON-CONVERGED fits (", if (DROP_NONCONVERGED) "EXCLUDED" else "kept", "):")
  print(nc, row.names = FALSE)
}
if (DROP_NONCONVERGED) fits <- fits[fits$converged, ]

write.csv(fits, file.path(OUT, "meta_group_estimates.csv"), row.names = FALSE)

# ---- random-effects meta-analysis, per treatment x effect -------------------
meta_one <- function(d) {
  # REML random-effects model; Knapp-Hartung adjustment for the CI on the mean
  m <- try(rma(yi = estimate, sei = se, data = d, method = "REML", test = "knha"),
           silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  ci_t2 <- try(confint(m), silent = TRUE)
  t2lo <- t2hi <- NA_real_
  if (!inherits(ci_t2, "try-error")) {
    t2lo <- ci_t2$random["tau^2", "ci.lb"]; t2hi <- ci_t2$random["tau^2", "ci.ub"]
  }
  data.frame(
    model = d$model[1], treatment = d$treatment[1], effect = d$effect[1],
    k = m$k,
    pooled = as.numeric(m$beta), pooled_se = m$se,
    ci_lb = m$ci.lb, ci_ub = m$ci.ub, p_pooled = m$pval,
    tau2 = m$tau2, tau2_lb = t2lo, tau2_ub = t2hi, tau = sqrt(m$tau2),
    I2 = m$I2, H2 = m$H2,
    Q = m$QE, Q_df = m$k - 1, Q_p = m$QEp,
    # tau relative to the pooled effect: I^2 says heterogeneity is DETECTABLE,
    # tau_rel says whether it is LARGE relative to the effect itself. An effect
    # can have I^2 ~ 80% with a negligible tau (precise estimates, tiny spread).
    tau_rel = sqrt(m$tau2) / abs(as.numeric(m$beta)),
    mean_se = mean(d$se), sd_estimate = sd(d$estimate),
    stringsAsFactors = FALSE
  )
}

key <- interaction(fits$model, fits$treatment, fits$effect, drop = TRUE)
res <- do.call(rbind, lapply(split(fits, key), meta_one))
res <- res[order(res$model, res$effect, res$treatment), ]

for (m in unique(res$model))
  write.csv(res[res$model == m, ], file.path(OUT, paste0("meta_", sub("_model", "", m), ".csv")),
            row.names = FALSE)
write.csv(res, file.path(OUT, "meta_all.csv"), row.names = FALSE)

# ---- report -----------------------------------------------------------------
cat("\n", strrep("=", 108), "\n",
    "RANDOM-EFFECTS META-ANALYSIS BY TREATMENT  (heterogeneity between groups)\n",
    strrep("=", 108), "\n", sep = "")
for (m in unique(res$model)) {
  cat("\n########", toupper(sub("_model", "", m)), "MODEL ########\n")
  for (tr in unique(res$treatment)) {
    d <- res[res$model == m & res$treatment == tr, ]
    if (!nrow(d)) next
    d <- d[order(-d$I2), ]
    cat("\n---", tr, "---\n")
    print(data.frame(
      effect = d$effect, k = d$k,
      pooled = sprintf("%8.4f", d$pooled),
      CI95 = sprintf("[%7.4f,%7.4f]", d$ci_lb, d$ci_ub),
      tau2 = sprintf("%.5f", d$tau2),
      tau = sprintf("%.4f", d$tau),
      I2 = sprintf("%5.1f%%", d$I2),
      Q_p = format.pval(d$Q_p, digits = 2, eps = 1e-4)
    ), row.names = FALSE)
  }
}

cat("\n", strrep("=", 108), "\nRANDOM-EFFECT CANDIDATES (ranked by mean I^2 across treatments)\n",
    strrep("=", 108), "\n", sep = "")
rank <- do.call(rbind, lapply(split(res, list(res$model, res$effect), drop = TRUE), function(d)
  data.frame(model = d$model[1], effect = d$effect[1],
             mean_I2 = mean(d$I2), min_I2 = min(d$I2), max_I2 = max(d$I2),
             mean_tau = mean(d$tau), mean_tau_rel = mean(d$tau_rel),
             n_treat_Qsig = sum(d$Q_p < 0.05), stringsAsFactors = FALSE)))
rank <- rank[order(rank$model, -rank$mean_I2), ]
write.csv(rank, file.path(OUT, "meta_candidates.csv"), row.names = FALSE)
print(data.frame(
  model = sub("_model", "", rank$model), effect = rank$effect,
  mean_I2 = sprintf("%5.1f%%", rank$mean_I2),
  range_I2 = sprintf("%.0f-%.0f%%", rank$min_I2, rank$max_I2),
  mean_tau = sprintf("%.4f", rank$mean_tau),
  tau_rel = sprintf("%6.2f", rank$mean_tau_rel),
  Q_sig = paste0(rank$n_treat_Qsig, "/3")
), row.names = FALSE)

# ---- forest plots -----------------------------------------------------------
PLOTS <- file.path(OUT, "plots"); dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)
for (m in unique(fits$model)) {
  for (ef in unique(fits$effect[fits$model == m])) {
    png(file.path(PLOTS, sprintf("forest_%s_%s.png", sub("_model", "", m), ef)),
        width = 2400, height = 2000, res = 200)
    op <- par(mfrow = c(3, 1), mar = c(4, 4, 2.5, 2))
    for (tr in unlist(TREATMENT_LAB)) {
      d <- fits[fits$model == m & fits$effect == ef & fits$treatment == tr, ]
      if (!nrow(d)) next
      mm <- try(rma(yi = estimate, sei = se, data = d, method = "REML", test = "knha"),
                silent = TRUE)
      if (inherits(mm, "try-error")) next
      forest(mm, slab = sub("^[^_]*_", "", d$group), header = c(tr, "estimate [95% CI]"),
             cex = 0.55, col = "#2a78d6", colout = "#22222f", efac = 0.4)
      mtext(sprintf("%s  |  I2 = %.1f%%   tau = %.4f   Q p = %s",
                    ef, mm$I2, sqrt(mm$tau2), format.pval(mm$QEp, digits = 2, eps = 1e-4)),
            side = 3, line = 0.4, cex = 0.62, col = "#6b6b76")
    }
    par(op); dev.off()
  }
}
message("wrote forest plots to ", PLOTS)

cat("\nInterpretation: I^2 is the share of variance between groups rather than within.",
    "\n  >75% = high heterogeneity -> strong candidate for a group-level random effect",
    "\n  tau  = between-group SD on the coefficient scale -> a prior scale for that random effect\n")
cat("\nwrote ", OUT, "/meta_all.csv, meta_choice.csv, meta_rate.csv, meta_group_estimates.csv\n", sep = "")
