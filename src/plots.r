#!/usr/bin/env Rscript
# =============================================================================
# Plot choice + rate results from the saved fits.
#   cd /home/milky/droso-pipe/6_dynam_models && Rscript plot_results.r
# Writes PNGs to res/choice/plots/ and res/rate/plots/.
# =============================================================================

suppressPackageStartupMessages({
  library(posterior); library(bayesplot); library(ggplot2)
})
theme_set(theme_minimal(base_size = 12))
bayesplot::color_scheme_set("brightblue")

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
res     <- file.path(PROJECT, "res")

save_png <- function(p, path, w = 9, h = 10) {
  ggplot2::ggsave(path, p, width = w, height = h, dpi = 300, limitsize = FALSE)
  message("wrote ", path)
}
q2 <- function(...) posterior::quantile2(..., probs = c(0.05, 0.95))  # -> q5, q95

# ===========================================================================
# CHOICE
# ===========================================================================
plot_choice <- function() {
  outdir <- file.path(res, "choice", "plots"); dir.create(outdir, FALSE, TRUE)
  fit   <- readRDS(file.path(res, "choice", "fit_choice.rds"))
  dsr   <- readRDS(file.path(res, "choice", "data_stan_re_choice.rds"))
  draws <- fit$draws()
  # X_choice was renamed to X_choice_raw in the pipeline; name it exactly so this
  # doesn't silently depend on `$` partial matching.
  xcols <- colnames(dsr$data_stan$X_choice_raw)
  # beta_choice covers every column EXCEPT the random-effect one (V1_choice),
  # whose coefficient the Stan model carries in mu/gamma instead.
  pnames <- xcols[-dsr$data_stan$V1_choice]
  C <- dsr$data_stan$C
  treatments <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")[seq_len(C)]
  message("choice: ", length(pnames), " fixed effects x ", C, " treatments, ",
          posterior::ndraws(draws), " draws")

  # --- fixed-effect intervals (beta_choice[c, p]: c = treatment, p = predictor) ---
  bs <- summarise_draws(subset_draws(draws, variable = "beta_choice"), median, q2)
  m  <- regmatches(bs$variable, regexec("beta_choice\\[(\\d+),(\\d+)\\]", bs$variable))
  ci <- as.integer(sapply(m, function(z) if (length(z) == 3) z[2] else NA))
  pj <- as.integer(sapply(m, function(z) if (length(z) == 3) z[3] else NA))
  bs <- bs[!is.na(ci), ]; ci <- ci[!is.na(ci)]; pj <- pj[!is.na(pj)]
  stopifnot(max(pj) == length(pnames))   # labels must line up with the draws
  bs$treatment <- factor(treatments[ci], levels = treatments)
  bs$label     <- factor(pnames[pj], levels = rev(pnames))
  p <- ggplot(bs, aes(median, label, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
    labs(x = "coefficient (90% CrI)", y = NULL, colour = NULL,
         title = "Choice model — fixed effects by treatment") +
    theme(legend.position = "top")
  save_png(p, file.path(outdir, "choice_beta_intervals.png"), 10,
           max(6, length(pnames) * 0.55))

  # --- per-group random effects (gamma) ---
  if ("gamma[1]" %in% variables(draws)) {
    gi <- readRDS(file.path(res, "choice", "groups_info.rds"))
    gs <- summarise_draws(subset_draws(draws, variable = "gamma"), median, q2)
    gs$label <- factor(gi$group[seq_len(nrow(gs))],
                       levels = gi$group[seq_len(nrow(gs))][order(gs$median)])
    pg <- ggplot(gs, aes(median, label)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0) + geom_point(size = 1.2) +
      labs(x = "random effect (90% CrI)", y = NULL, title = "Choice model — group random effects")
    save_png(pg, file.path(outdir, "choice_gamma_intervals.png"), 8, max(6, nrow(gs) * 0.18))
  }

  # --- diagnostics ---
  sm <- fit$summary()
  nd <- posterior::ndraws(draws)
  save_png(mcmc_rhat(setNames(sm$rhat, sm$variable)) + ggtitle("Choice — Rhat"),
           file.path(outdir, "choice_rhat.png"), 7, 10)
  save_png(mcmc_neff(setNames(sm$ess_bulk / nd, sm$variable)) + ggtitle("Choice — ESS ratio"),
           file.path(outdir, "choice_ess.png"), 7, 10)
  save_png(mcmc_trace(draws, regex_pars = "beta_choice\\[1,[1-6]\\]") + ggtitle("Choice — traces"),
           file.path(outdir, "choice_trace.png"), 9, 8)
  message("choice: worst Rhat = ", round(max(sm$rhat, na.rm = TRUE), 4),
          " | min ESS ratio = ", round(min(sm$ess_bulk / nd, na.rm = TRUE), 3))
}

# ===========================================================================
# RATE   (beta_rate[c, p]: c = treatment category, p = predictor)
# ===========================================================================
plot_rate <- function() {
  outdir <- file.path(res, "rate", "plots"); dir.create(outdir, FALSE, TRUE)
  fit   <- readRDS(file.path(res, "rate", "fit_rate.rds"))
  dsr   <- readRDS(file.path(res, "rate", "data_stan_re_rate.rds"))
  draws <- fit$draws()
  pnames <- colnames(dsr$data_stan$X_rate_raw)
  C <- dsr$data_stan$C
  treatments <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")[seq_len(C)]
  message("rate: ", length(pnames), " predictors x ", C, " treatments, ",
          posterior::ndraws(draws), " draws")

  # --- beta_rate by treatment (original scale, back-transformed) ---
  s <- summarise_draws(subset_draws(draws, variable = "beta_rate"), median, q2)
  m  <- regmatches(s$variable, regexec("beta_rate\\[(\\d+),(\\d+)\\]", s$variable))
  ci <- as.integer(sapply(m, function(z) if (length(z) == 3) z[2] else NA))
  pj <- as.integer(sapply(m, function(z) if (length(z) == 3) z[3] else NA))
  s  <- s[!is.na(ci), ]; ci <- ci[!is.na(ci)]; pj <- pj[!is.na(pj)]
  s$treatment <- factor(treatments[ci], levels = treatments)
  s$predictor <- factor(pnames[pj], levels = rev(pnames))
  pd <- position_dodge(width = 0.6)
  p <- ggplot(s, aes(median, predictor, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0, position = pd) +
    geom_point(size = 1.8, position = pd) +
    labs(x = "rate coefficient (90% CrI)", y = NULL, colour = NULL,
         title = "Rate model — effects by treatment") +
    theme(legend.position = "top")
  save_png(p, file.path(outdir, "rate_beta_by_treatment.png"), 10, max(6, length(pnames) * 0.55))

  # --- treatment contrasts vs CS_10D/old baseline (computed from draws) ---
  dd <- as_draws_df(subset_draws(draws, variable = "beta_rate"))
  rows <- do.call(rbind, lapply(seq_along(pnames), function(p) {
    base <- dd[[sprintf("beta_rate[1,%d]", p)]]
    out <- list()
    if (C >= 2) out[[length(out) + 1]] <- data.frame(
      predictor = pnames[p], contrast = "isolated \u2212 old",
      median = median(dd[[sprintf("beta_rate[2,%d]", p)]] - base),
      q5 = quantile(dd[[sprintf("beta_rate[2,%d]", p)]] - base, .05),
      q95 = quantile(dd[[sprintf("beta_rate[2,%d]", p)]] - base, .95))
    if (C >= 3) out[[length(out) + 1]] <- data.frame(
      predictor = pnames[p], contrast = "young \u2212 old",
      median = median(dd[[sprintf("beta_rate[3,%d]", p)]] - base),
      q5 = quantile(dd[[sprintf("beta_rate[3,%d]", p)]] - base, .05),
      q95 = quantile(dd[[sprintf("beta_rate[3,%d]", p)]] - base, .95))
    do.call(rbind, out)
  }))
  rows$predictor <- factor(rows$predictor, levels = rev(pnames))
  pc <- ggplot(rows, aes(median, predictor, colour = contrast)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0, position = position_dodge(0.6)) +
    geom_point(size = 1.8, position = position_dodge(0.6)) +
    labs(x = "difference in rate coefficient (90% CrI)", y = NULL, colour = NULL,
         title = "Rate model — treatment contrasts vs CS_10D (old)") +
    theme(legend.position = "top")
  save_png(pc, file.path(outdir, "rate_treatment_contrasts.png"), 10, max(6, length(pnames) * 0.55))

  # --- diagnostics ---
  sm <- fit$summary()
  nd <- posterior::ndraws(draws)
  save_png(mcmc_rhat(setNames(sm$rhat, sm$variable)) + ggtitle("Rate — Rhat"),
           file.path(outdir, "rate_rhat.png"), 7, 10)
  save_png(mcmc_neff(setNames(sm$ess_bulk / nd, sm$variable)) + ggtitle("Rate — ESS ratio"),
           file.path(outdir, "rate_ess.png"), 7, 10)
  save_png(mcmc_trace(draws, regex_pars = "beta_rate\\[1,[1-6]\\]|sigma") + ggtitle("Rate — traces"),
           file.path(outdir, "rate_trace.png"), 9, 8)
  message("rate: worst Rhat = ", round(max(sm$rhat, na.rm = TRUE), 4),
          " | min ESS ratio = ", round(min(sm$ess_bulk / nd, na.rm = TRUE), 3))
}

# ---- run both, don't let one failure stop the other ------------------------
tryCatch(plot_choice(), error = function(e) message("CHOICE plotting failed: ", conditionMessage(e)))
tryCatch(plot_rate(),   error = function(e) message("RATE plotting failed: ",   conditionMessage(e)))
message("\nDone. PNGs in res/choice/plots/ and res/rate/plots/")