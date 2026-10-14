#!/usr/bin/env Rscript
# =============================================================================
# Aggregate the posterior predictive check with credible intervals.
#   Rscript src/5f_aggregate.r
# Reads res/gof/gof_raw.rds (written by 5e_gof_run.r) and writes:
#   gof_by_treatment.csv  treatment-level: observed vs predictive 90% + 95% CrI
#   gof_by_group.csv      per-group detail within each treatment
#
# NOTE ON TERMINOLOGY: the intervals here are posterior PREDICTIVE intervals for
# a replicated statistic, not credible intervals on a model parameter. Parameter
# CrIs live in res/posterior/*_summary.csv (q2.5/q5/q95/q97.5).
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
OUT <- file.path("res", "gof")
TREATMENTS <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")

results <- readRDS(file.path(OUT, "gof_raw.rds"))
cat("groups loaded:", length(results), "\n")
stat_names <- names(results[[1]]$obs)
n_draws <- length(results[[1]]$sim)
cat("draws per group:", n_draws, " statistics:", length(stat_names), "\n")

qq <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

# ---- per group --------------------------------------------------------------
by_group <- do.call(rbind, lapply(results, function(r) {
  do.call(rbind, lapply(stat_names, function(st) {
    v <- vapply(r$sim, function(s) as.numeric(s[[st]]), 0)
    o <- as.numeric(r$obs[[st]])
    data.frame(
      treatment = TREATMENTS[r$cc], group = r$group, statistic = st,
      observed = o,
      sim_median = median(v, na.rm = TRUE), sim_mean = mean(v, na.rm = TRUE),
      sim_q2.5 = qq(v, .025), sim_q5 = qq(v, .05),
      sim_q95 = qq(v, .95), sim_q97.5 = qq(v, .975),
      ppp = mean(v >= o, na.rm = TRUE),
      inside_90 = o >= qq(v, .05)  & o <= qq(v, .95),
      inside_95 = o >= qq(v, .025) & o <= qq(v, .975),
      n_draws = sum(!is.na(v)),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(by_group, file.path(OUT, "gof_by_group.csv"), row.names = FALSE)

# ---- per treatment ----------------------------------------------------------
# For each draw, average the statistic across the groups in the treatment, then
# summarise that distribution of treatment-level means.
by_treatment <- do.call(rbind, lapply(TREATMENTS, function(tr) {
  rs <- Filter(function(r) TREATMENTS[r$cc] == tr, results)
  if (!length(rs)) return(NULL)
  do.call(rbind, lapply(stat_names, function(st) {
    M <- vapply(rs, function(r) vapply(r$sim, function(s) as.numeric(s[[st]]), 0),
                numeric(n_draws))                      # draws x groups
    per_draw <- rowMeans(M, na.rm = TRUE)
    o <- mean(vapply(rs, function(r) as.numeric(r$obs[[st]]), 0), na.rm = TRUE)
    data.frame(
      treatment = tr, statistic = st, n_groups = length(rs),
      observed = o,
      sim_median = median(per_draw, na.rm = TRUE), sim_mean = mean(per_draw, na.rm = TRUE),
      sim_q2.5 = qq(per_draw, .025), sim_q5 = qq(per_draw, .05),
      sim_q95 = qq(per_draw, .95), sim_q97.5 = qq(per_draw, .975),
      ppp = mean(per_draw >= o, na.rm = TRUE),
      inside_90 = o >= qq(per_draw, .05)  & o <= qq(per_draw, .95),
      inside_95 = o >= qq(per_draw, .025) & o <= qq(per_draw, .975),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(by_treatment, file.path(OUT, "gof_by_treatment.csv"), row.names = FALSE)

cat("\n", strrep("=", 100), "\nPOSTERIOR PREDICTIVE CHECK BY TREATMENT (95% predictive interval)\n",
    strrep("=", 100), "\n", sep = "")
for (tr in TREATMENTS) {
  d <- by_treatment[by_treatment$treatment == tr, ]
  if (!nrow(d)) next
  cat("\n---", tr, "(", d$n_groups[1], "groups ) ---\n")
  out <- data.frame(
    statistic = d$statistic,
    observed  = signif(d$observed, 4),
    predicted = signif(d$sim_median, 4),
    CrI95     = sprintf("[%.4g, %.4g]", d$sim_q2.5, d$sim_q97.5),
    ppp       = round(d$ppp, 3),
    flag      = ifelse(d$inside_95, "", "  <-- OUTSIDE")
  )
  print(out, row.names = FALSE)
}
cat("\ncoverage: ", sum(by_treatment$inside_95), "/", nrow(by_treatment),
    " treatment-level statistics inside their 95% predictive interval\n", sep = "")
cat("wrote gof_by_treatment.csv and gof_by_group.csv to ", OUT, "\n", sep = "")
