#!/usr/bin/env Rscript
# =============================================================================
# Forward simulation of the fitted Bayesian DyNAM, using the VALIDATED
# statistic engine in 5a_stat_engine.r (proven to reproduce goldfish exactly:
# choice 26/26 columns, rate 9/9 columns).
#
# Per posterior draw s, group a in treatment c:
#   log rate_i(t) = alpha[c,s] + gamma_rate[a,s] + x_rate_i(t) . beta_rate[c,,s]
#   utility_ij(t) = x_choice_ij(t) . c(beta_choice[c,,s], mu[c,s] + gamma_ch[a,s])
#
# Rates are piecewise constant between covariate changes (integer times) and
# network changes (events), so events are drawn exactly: sample a waiting time
# from Exp(sum of rates); if it lands before the next covariate boundary an
# event occurs, otherwise advance to the boundary and recompute.
#
# Covariates are REPLAYED from the observed trajectories (matching how the
# model was fitted); simulation runs over the observed time horizon.
# =============================================================================

source(file.path(dirname(sys.frame(1)$ofile %||% "src/x"), "5a_stat_engine.r"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- assemble everything one group needs ------------------------------------
build_group <- function(a, gi, base = ".") {
  grp <- as.character(gi$group[a]); trt <- as.character(gi$treatment[a])
  ev <- read.csv(file.path(base, "data/edgelists", trt, paste0(grp, ".csv")))
  ev <- ev[order(ev$time), c("time", "sender", "receiver", "increment")]
  labels <- as.character(sort(unique(c(ev$sender, ev$receiver))))
  gcp <- file.path(base, "data/covariances", trt, grp)

  cov_files <- c(
    activity                               = "out_degree.csv",
    positive_influence                     = "positiveinfluence.csv",
    negative_influence                     = "negativeinfluence.csv",
    positive_inf_weighted                  = "positiveinf_weighted.csv",
    negative_inf_weighted                  = "negativeinf_weighted.csv",
    distance_traveled_between_interactions = "distance_traveled_between_interactions.csv",
    number_of_flies_in_soc_space           = "number_of_flies_in_soc_space.csv",
    unique_partners_met_social_space       = "unique_partners_met_social_space.csv"
  )
  covs <- lapply(names(cov_files), function(nm) {
    df <- read_and_clean(cov_files[[nm]], gcp)
    # 3_rate_fit.r log1p-transformed distance AFTER gathering -> match it
    if (nm == "distance_traveled_between_interactions") df$replace <- log1p(df$replace)
    df
  })
  names(covs) <- names(cov_files)

  list(index = a, group = grp, treatment = trt, labels = labels,
       n = length(labels), obs = ev, covs = covs,
       T_max = max(ev$time), n_obs = nrow(ev))
}

# ---- one forward simulation --------------------------------------------------
simulate_group <- function(G, p_rate, p_choice,
                           spec_rate = rate_spec(), spec_choice = choice_spec(),
                           max_events = 50000) {
  n <- G$n; T_max <- G$T_max
  cs <- lapply(G$covs, new_cov_state, labels = G$labels)
  WS <- new_window_state(CHOICE_WINDOWS, n)

  # simulated event log, grown in blocks (the engine reads it to expire windows)
  cap <- 4096L
  st <- list(time = numeric(cap), s = integer(cap), r = integer(cap),
             inc = numeric(cap), n = 0L)
  push <- function(st, t, i, j) {
    if (st$n == length(st$time)) {
      st$time <- c(st$time, numeric(length(st$time)))
      st$s <- c(st$s, integer(length(st$s)))
      st$r <- c(st$r, integer(length(st$r)))
      st$inc <- c(st$inc, numeric(length(st$inc)))
    }
    k <- st$n + 1L
    st$time[k] <- t; st$s[k] <- i; st$r[k] <- j; st$inc[k] <- 1
    st$n <- k
    st
  }

  t <- 0
  while (t < T_max && st$n < max_events) {
    for (k in seq_along(WS)) WS[[k]] <- expire(WS[[k]], st, t)
    for (nm in names(cs))    cs[[nm]] <- advance_cov(cs[[nm]], t)
    cov_cur <- lapply(cs, `[[`, "cur")

    lambda <- exp(p_rate[1] + rate_stats(WS, cov_cur, spec_rate, n) %*% p_rate[-1])
    Lam <- sum(lambda)
    if (!is.finite(Lam) || Lam <= 0) break

    dt <- rexp(1, Lam)
    boundary <- floor(t) + 1                       # next covariate change
    if (t + dt >= boundary) { t <- boundary; next } # no event in this segment
    t <- t + dt

    i   <- sample.int(n, 1L, prob = as.numeric(lambda))
    alt <- seq_len(n)[-i]
    u   <- choice_stats(i, alt, WS, cov_cur, spec_choice) %*% p_choice
    u   <- u - max(u)
    j   <- alt[sample.int(length(alt), 1L, prob = exp(as.numeric(u)))]

    st <- push(st, t, i, j)
    for (k in seq_along(WS)) WS[[k]] <- add_event(WS[[k]], i, j, 1)
  }

  data.frame(time = st$time[seq_len(st$n)],
             sender   = G$labels[st$s[seq_len(st$n)]],
             receiver = G$labels[st$r[seq_len(st$n)]],
             stringsAsFactors = FALSE)
}

# ---- GoF statistics (same five as the original frequentist script) ----------
network_stats <- function(ev, labels) {
  N <- length(labels)
  if (nrow(ev) == 0)
    return(list(indeg = rep(0, N), outdeg = rep(0, N), tie_weights = numeric(0),
                recip_ratio = 0, w_density = 0, n_events = 0))
  indeg  <- as.numeric(table(factor(ev$receiver, levels = labels)))
  outdeg <- as.numeric(table(factor(ev$sender,   levels = labels)))
  tw <- aggregate(list(weight = rep(1, nrow(ev))),
                  by = list(sender = ev$sender, receiver = ev$receiver), FUN = sum)
  dyads <- paste(tw$sender, tw$receiver, sep = "->")
  rev_d <- paste(tw$receiver, tw$sender, sep = "->")
  list(indeg = indeg, outdeg = outdeg, tie_weights = tw$weight,
       recip_ratio = mean(rev_d %in% dyads),
       w_density = nrow(ev) / (N * (N - 1)),
       n_events = nrow(ev))
}

# ---- posterior parameter draws for a group ----------------------------------
make_params <- function(a, cc, s, DR) {
  list(
    rate = c(DR$rate[s, sprintf("alpha[%d]", cc)] + DR$rate[s, sprintf("gamma[%d]", a)],
             vapply(1:9, function(p) DR$rate[s, sprintf("beta_rate[%d,%d]", cc, p)], 0)),
    choice = c(vapply(1:25, function(p) DR$choice[s, sprintf("beta_choice[%d,%d]", cc, p)], 0),
               DR$choice[s, sprintf("mu[%d]", cc)] + DR$choice[s, sprintf("gamma[%d]", a)])
  )
}
