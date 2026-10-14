#!/usr/bin/env Rscript
# =============================================================================
# Effect-statistic engine for DyNAM forward simulation (pure R, incremental).
#
# This file defines the statistics ONLY. It is validated in 5b_validate_engine.r
# by replaying the observed events and checking it reproduces goldfish's own
# X_choice_raw / X_rate_raw column-for-column. Nothing simulates until it does.
#
# Design: for each window W keep a weighted adjacency matrix A_W holding the
# summed increments of events currently inside the window, plus a head pointer
# into the event list. Advancing time expires old events by subtracting them,
# so no statistic is ever recomputed by rescanning history.
# =============================================================================

# ---- windowed network state -------------------------------------------------

new_window_state <- function(windows, n) {
  lapply(windows, function(w) list(
    w    = w,
    head = 1L,                       # oldest event still inside the window
    A    = matrix(0, n, n)           # summed increments inside the window
  ))
}

# expire events that have dropped out of the window as time advances to `now`
expire <- function(st, events, now) {
  # strict: an event at exactly `now - w` is still INSIDE the window (matches goldfish)
  while (st$head <= events$n && events$time[st$head] < now - st$w) {
    h <- st$head
    st$A[events$s[h], events$r[h]] <- st$A[events$s[h], events$r[h]] - events$inc[h]
    st$head <- h + 1L
  }
  st
}

add_event <- function(st, s, r, inc) {
  st$A[s, r] <- st$A[s, r] + inc
  st
}

# ---- time-varying node covariates -------------------------------------------
# cov_df: data.frame(time, node, replace) as produced by read_and_clean().
# Returns a stepper holding current values, advanced monotonically with time.

new_cov_state <- function(cov_df, labels, init = 0) {
  ord <- order(cov_df$time)
  list(
    time = cov_df$time[ord],
    node = match(as.character(cov_df$node[ord]), labels),
    val  = cov_df$replace[ord],
    i    = 1L,
    cur  = setNames(rep(init, length(labels)), labels),
    n    = nrow(cov_df)
  )
}

advance_cov <- function(cs, now) {
  # strict: attribute changes take effect only for events STRICTLY after them
  while (cs$i <= cs$n && cs$time[cs$i] < now) {
    if (!is.na(cs$node[cs$i])) cs$cur[cs$node[cs$i]] <- cs$val[cs$i]
    cs$i <- cs$i + 1L
  }
  cs
}

# ---- statistics -------------------------------------------------------------
# All take the CURRENT (pre-event) state and return values for the alternatives.

# binary view of a windowed network
bin <- function(A) (A > 0) * 1

# choice statistics for sender i over alternatives alt (vector of node indices)
choice_stats <- function(i, alt, WS, cov_cur, spec) {
  out <- matrix(0, length(alt), nrow(spec))
  Bc <- vector("list", length(WS)); Ac <- vector("list", length(WS))
  for (k in seq_along(WS)) { Ac[[k]] <- WS[[k]]$A; Bc[[k]] <- bin(Ac[[k]]) }

  for (col in seq_len(nrow(spec))) {
    k    <- spec$win_ix[col]
    kind <- spec$kind[col]
    out[, col] <- switch(kind,
      indeg      = colSums(Bc[[k]])[alt],
      outdeg     = rowSums(Bc[[k]])[alt],
      trans      = (Bc[[k]] %*% Bc[[k]])[i, alt],
      inertia_w  = Ac[[k]][i, alt],
      inertia_u  = Bc[[k]][i, alt],
      recip      = Bc[[k]][alt, i],
      alter      = cov_cur[[spec$attr[col]]][alt],
      sim        = -abs(cov_cur[[spec$attr[col]]][i] - cov_cur[[spec$attr[col]]][alt]),
      stop("unknown kind: ", kind)
    )
  }
  out
}

# rate statistics for every actor (ego-type)
rate_stats <- function(WS, cov_cur, spec, n) {
  out <- matrix(0, n, nrow(spec))
  for (col in seq_len(nrow(spec))) {
    k <- spec$win_ix[col]
    out[, col] <- switch(spec$kind[col],
      indeg  = colSums(bin(WS[[k]]$A)),
      outdeg = rowSums(bin(WS[[k]]$A)),
      ego    = cov_cur[[spec$attr[col]]],
      stop("unknown kind: ", spec$kind[col])
    )
  }
  out
}

# ---- effect specifications --------------------------------------------------
# Term order MUST match the fitted design matrices exactly; parameter vectors
# are positional. Choice: 25 fitted terms, then the random-effect term last.

CHOICE_WINDOWS <- c(96, 288, 864, 2592)

choice_spec <- function() {
  s <- function(kind, win, attr = NA_character_)
    data.frame(kind = kind, win = win, attr = attr, stringsAsFactors = FALSE)
  spec <- rbind(
    do.call(rbind, lapply(CHOICE_WINDOWS, function(w) s("indeg",  w))),
    do.call(rbind, lapply(CHOICE_WINDOWS, function(w) s("outdeg", w))),
    do.call(rbind, lapply(CHOICE_WINDOWS, function(w) s("trans",  w))),
    do.call(rbind, lapply(CHOICE_WINDOWS, function(w) s("inertia_w", w))),
    do.call(rbind, lapply(CHOICE_WINDOWS, function(w) s("inertia_u", w))),   # all 4 fixed (288 no longer the RE)
    do.call(rbind, lapply(c(96, 864, 2592), function(w) s("recip",  w))),    # recip 288 is now the RE
    s("alter", NA, "number_of_flies_in_soc_space"),
    s("sim",   NA, "activity"),
    s("recip", 288)                           # <- random-effect term, LAST (recip window 288)
  )
  spec$win_ix <- match(spec$win, CHOICE_WINDOWS)
  spec
}

rate_spec <- function() {
  s <- function(kind, win, attr = NA_character_)
    data.frame(kind = kind, win = win, attr = attr, stringsAsFactors = FALSE)
  spec <- rbind(
    s("indeg", 96), s("outdeg", 96),
    s("ego", NA, "positive_influence"),
    s("ego", NA, "negative_influence"),
    s("ego", NA, "positive_inf_weighted"),
    s("ego", NA, "negative_inf_weighted"),
    s("ego", NA, "distance_traveled_between_interactions"),
    s("ego", NA, "number_of_flies_in_soc_space"),
    s("ego", NA, "unique_partners_met_social_space")
  )
  spec$win_ix <- match(spec$win, CHOICE_WINDOWS)
  spec
}
