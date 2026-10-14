#!/usr/bin/env Rscript
# =============================================================================
# Posterior predictive check figures, faceted by treatment.
#   Rscript src/5g_gof_plots.r
# Writes PNGs to res/gof/plots/.
#
# Colour: two series only -- model (predictive) vs observed -- assigned from the
# reference categorical palette slots 1 and 2 in fixed order.
# =============================================================================

PROJECT <- "/home/milky/droso-pipe/6_dynam_models"
setwd(PROJECT)
suppressPackageStartupMessages({ library(ggplot2) })

OUT <- file.path("res", "gof", "plots"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TREATMENTS <- c("CS_10D (old)", "Cs_5DIZ (isolated)", "CsCh (young)")

SIM_COL <- "#2a78d6"   # categorical slot 1 - model / predictive
OBS_COL <- "#eb6834"   # categorical slot 2 - observed
INK     <- "#22222f"; MUTED <- "#6b6b76"

theme_gof <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
    axis.title = element_text(colour = MUTED, size = 9),
    axis.text  = element_text(colour = INK),
    strip.text = element_text(colour = INK, face = "bold", size = 9.5),
    plot.title = element_text(colour = INK, face = "bold", size = 13),
    plot.subtitle = element_text(colour = MUTED, size = 9.5),
    plot.caption  = element_text(colour = MUTED, size = 8, hjust = 0),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(colour = INK, size = 9)
  )
save_png <- function(p, f, w, h) {
  ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "#fcfcfb")
  message("wrote ", file.path(OUT, f))
}

results <- readRDS(file.path("res", "gof", "gof_raw.rds"))
message("groups: ", length(results), " | draws: ", length(results[[1]]$sim))

# =============================================================================
# FIGURE 1 - scalar statistics: observed vs predictive interval
# =============================================================================
tb <- read.csv(file.path("res", "gof", "gof_by_treatment.csv"), stringsAsFactors = FALSE)

# drop statistics that are saturated (degenerate interval in observed AND model):
# they carry no information and would read as spurious misfit
sat <- tapply(tb$sim_q97.5 - tb$sim_q2.5, tb$statistic, function(x) all(x < 1e-6))
drop <- names(sat)[sat]
if (length(drop)) message("dropping saturated statistics: ", paste(drop, collapse = ", "))
tb <- tb[!tb$statistic %in% drop, ]

lbl <- c(n_events = "events", w_density = "weighted density", tie_w_mean = "tie weight (mean)",
         tie_w_sd = "tie weight (SD)", tie_w_max = "tie weight (max)",
         indeg_sd = "in-degree SD", outdeg_sd = "out-degree SD",
         indeg_gini = "in-degree Gini", outdeg_gini = "out-degree Gini",
         iet_median = "inter-event time (median)", iet_cv = "inter-event time (CV)",
         burstiness = "burstiness", repeat_frac = "repeated dyad fraction",
         two_paths = "two-paths", recip_ratio = "reciprocity", transitivity = "transitivity")
tb$stat_lab <- ifelse(is.na(lbl[tb$statistic]), tb$statistic, lbl[tb$statistic])
tb$treatment <- factor(tb$treatment, levels = TREATMENTS)
tb$stat_lab <- factor(tb$stat_lab, levels = unique(tb$stat_lab))

p1 <- ggplot(tb, aes(y = treatment)) +
  geom_linerange(aes(xmin = sim_q2.5, xmax = sim_q97.5, colour = "model (95% predictive)"),
                 linewidth = 0.9) +
  geom_linerange(aes(xmin = sim_q5, xmax = sim_q95, colour = "model (95% predictive)"),
                 linewidth = 2.2, alpha = 0.45) +
  geom_point(aes(x = sim_median, colour = "model (95% predictive)"), size = 2.4) +
  geom_point(aes(x = observed, colour = "observed"), size = 3.1, shape = 18) +
  facet_wrap(~ stat_lab, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("model (95% predictive)" = SIM_COL, "observed" = OBS_COL)) +
  scale_y_discrete(limits = rev) +
  labs(title = "Posterior predictive check by treatment",
       subtitle = "Diamond = observed. Line = model's 95% predictive interval for the treatment-level mean (thick band 90%).",
       x = NULL, y = NULL,
       caption = paste0("100 posterior draws x 20 groups per treatment. Forward simulation over the observed time horizon.\n",
                        if (length(drop)) paste0("Saturated statistics omitted: ", paste(drop, collapse = ", "), ".") else "")) +
  theme_gof
save_png(p1, "gof_scalar_intervals.png", 10, 11)

# =============================================================================
# FIGURE 2 - distributions pooled within treatment
# =============================================================================
# For each draw: pool the quantity across the groups of a treatment, bin it,
# and take per-bin quantiles ACROSS draws -> predictive ribbon.
# Observed event times are integers, so observed inter-event times are integer-
# valued while the continuous-time model produces fractional ones. Comparing them
# on arbitrary bins produces a spurious comb. `snap` puts simulated times back on
# the integer grid the model was fitted on, making the comparison like-for-like.
snap_iet <- function(x) { d <- diff(round(cumsum(x))); d[d > 0] }

ribbon_for <- function(field, breaks, treat, snap = FALSE) {
  rs <- Filter(function(r) TREATMENTS[r$cc] == treat, results)
  nd <- length(rs[[1]]$sim_dist)
  mids <- head(breaks, -1) + diff(breaks) / 2
  bw   <- diff(breaks)          # divide by bin width -> density, so a change in
                                # binning does not create a spurious step
  M <- vapply(seq_len(nd), function(s) {
    v <- unlist(lapply(rs, function(r) {
      z <- r$sim_dist[[s]][[field]]; if (snap) snap_iet(z) else z
    }))
    as.numeric(table(cut(v, breaks, include.lowest = TRUE))) / bw
  }, numeric(length(mids)))
  obs <- unlist(lapply(rs, function(r) r$obs_dist[[field]]))
  data.frame(
    treatment = treat, x = mids,
    lo95 = apply(M, 1, quantile, .025), lo90 = apply(M, 1, quantile, .05),
    med  = apply(M, 1, median),
    hi90 = apply(M, 1, quantile, .95),  hi95 = apply(M, 1, quantile, .975),
    observed = as.numeric(table(cut(obs, breaks, include.lowest = TRUE))) / bw
  )
}

dist_panel <- function(field, breaks, title, xlab, logx = FALSE, file, h = 3.6,
                       snap = FALSE, sub = NULL) {
  d <- do.call(rbind, lapply(TREATMENTS, function(t) ribbon_for(field, breaks, t, snap)))
  # trim trailing empty bins per treatment so panels are not mostly whitespace
  d <- do.call(rbind, lapply(split(d, d$treatment), function(z) {
    nz <- which(z$observed > 0 | z$hi95 > 0)
    if (!length(nz)) return(z)
    z[seq(max(1, min(nz) - 1), min(nrow(z), max(nz) + 1)), ]
  }))
  d$treatment <- factor(d$treatment, levels = TREATMENTS)
  p <- ggplot(d, aes(x = x)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = "model (95% predictive)"), alpha = 0.22) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90, fill = "model (95% predictive)"), alpha = 0.28) +
    geom_line(aes(y = med, colour = "model (95% predictive)"), linewidth = 0.7) +
    geom_line(aes(y = observed, colour = "observed"), linewidth = 0.8) +
    geom_point(aes(y = observed, colour = "observed"), size = 1.5) +
    facet_wrap(~ treatment, nrow = 1, scales = "free") +
    scale_colour_manual(values = c("model (95% predictive)" = SIM_COL, "observed" = OBS_COL)) +
    scale_fill_manual(values  = c("model (95% predictive)" = SIM_COL), guide = "none") +
    labs(title = title, x = xlab, y = "density (count / bin width)",
         subtitle = sub %||% "Pooled across the 20 groups of each treatment; ribbon = predictive interval across 100 draws") +
    theme_gof
  if (logx) p <- p + scale_x_log10()
  save_png(p, file, 11, h)
}

all_field <- function(f) unlist(lapply(results, function(r) r$obs_dist[[f]]))

deg_max <- max(all_field("indeg"), all_field("outdeg"))
dist_panel("indeg",  seq(0, deg_max * 1.4, length.out = 26),
           "In-degree distribution", "in-degree (events received)", FALSE, "gof_dist_indegree.png")
dist_panel("outdeg", seq(0, deg_max * 1.4, length.out = 26),
           "Out-degree distribution", "out-degree (events sent)", FALSE, "gof_dist_outdegree.png")

tw_max <- max(all_field("tie_weights"))
dist_panel("tie_weights", seq(0, tw_max * 1.4, length.out = 26),
           "Tie-weight distribution", "events per directed dyad", FALSE, "gof_dist_tieweight.png")

# integer-aligned bins: each of 1..20 gets its own bin, then log-spaced above
iet <- all_field("iet"); iet <- iet[iet > 0]
iet_breaks <- c(seq(0.5, 20.5, by = 1),
                10^seq(log10(21), log10(max(iet) * 1.5), length.out = 12))
dist_panel("iet", iet_breaks, "Inter-event time distribution",
           "inter-event time (log scale)", TRUE, "gof_dist_iet.png", snap = TRUE,
           sub = paste("Pooled across the 20 groups per treatment; ribbon = predictive interval across 100 draws.",
                       "Simulated times snapped to the integer grid so both are integer-valued."))

message("\nDone. PNGs in ", OUT)
