# =============================================================================
# DyNAM choice model pipeline — Drosophila social interaction networks
#
# Targets:
#   - goldfish @ refactor/rate_prep  (native snake_case API: make_data,
#     gather_model_data, estimate_dynam — fixed effects, no random-effects layer)
#
# Install (run once, then comment out):
#   install.packages(c("here","remotes","glue","foreach","doParallel",
#                      "broom","corrplot"))
#   remotes::install_github("stocnet/goldfish@refactor/rate_prep")
#
# Pipeline stages:
#   0. Repeat-count diagnostic     -> decide weighted vs unweighted empirically
#   1. Per-group preprocessing     -> choice design matrices (.rds)
#   2. Collinearity diagnostics    -> high-correlation plot
#   3. Native DyNAM choice fits    -> per-group coefficient CSVs (+ weighted refit)
# =============================================================================

suppressPackageStartupMessages({
  library(goldfish)
  library(here)
  library(glue)
  library(foreach)
  library(doParallel)
})

# -----------------------------------------------------------------------------
# read_and_clean(): load a covariate change-event CSV.
# Defined inline (previously lived in src/utils.R, which is not present). The
# covariate CSVs are already in goldfish change-event format — columns
# time, actor, value — so this is a thin, validated read. Kept as a named
# function so parallel workers can call it without sourcing an external file.
# -----------------------------------------------------------------------------
read_and_clean <- function(filename, directory) {
  fpath <- file.path(directory, filename)
  df <- read.csv(fpath, stringsAsFactors = FALSE)

  # goldfish node change-events require columns: time, node, replace.
  # Our CSVs ship as time, actor, value -> rename to what link_events expects.
  if ("actor" %in% names(df)) names(df)[names(df) == "actor"] <- "node"
  if ("value" %in% names(df)) names(df)[names(df) == "value"] <- "replace"

  needed <- c("time", "node", "replace")
  if (!all(needed %in% names(df))) {
    stop(glue(
      "unexpected columns in {filename}: got [{paste(names(df), collapse = ', ')}], ",
      "need [time, node, replace] (after renaming actor->node, value->replace)"
    ), call. = FALSE)
  }
  df <- df[, needed]
  df$time    <- as.numeric(df$time)
  df$replace <- as.numeric(df$replace)

  # Apply a read-time transform if this file is listed in cov_transform
  # (e.g. log1p for the heavily right-skewed distance covariate).
  if (exists("cov_transform") && filename %in% names(cov_transform)) {
    tf <- cov_transform[[filename]]
    df$replace <- switch(tf,
      log1p = log1p(df$replace),
      log   = log(df$replace),
      df$replace  # unknown transform -> leave as-is
    )
  }

  df <- df[order(df$time), ]
  df
}

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
cfg <- list(
  treatments  = c("CS_10D", "Cs_5DIZ", "CsCh"),
  base_path   = "/home/milky/droso-pipe/6_dynam_models",
  cor_thresh  = 0.75,
  # Cores for the parallel workers. On a SHARED node, requesting too many
  # triggers "fork: Resource temporarily unavailable" (you hit the per-user
  # process limit; check `ulimit -u`). Each worker can also spawn threads, so
  # stay well under the machine total. Capped at 12 by default; raise only if
  # the node is idle and your process limit is high.
  n_cores     = min(12L, max(1L, parallel::detectCores() - 2L)),
  # set to an integer to process only the first N groups (debugging); NULL = all
  limit_groups = NULL
)

paths <- within(list(), {
  edgelists      <- file.path(cfg$base_path, "data", "edgelists")
  covariances    <- file.path(cfg$base_path, "data", "covariances")
  soc_space      <- file.path(cfg$base_path, "data", "soc_space_matrix")
  results        <- file.path(cfg$base_path, "res")
  tmp            <- file.path(cfg$base_path, "tmp")
})
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tmp,     showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(paths$results, "parallel_execution_log.txt")

# Treatment -> which age/isolation indicator flag is switched on for that group.
treatment_flag <- c(CsCh = "young", CS_10D = "old", Cs_5DIZ = "isolated")

# Covariance files linked onto node attributes: file name -> attribute name.
cov_attr_map <- c(
  in_degree.csv                              = "popularity",              # sim(popularity)
  positiveinfluence.csv                      = "positive_influence",      # alter/ego
  negativeinfluence.csv                      = "negative_influence",      # alter/ego
  number_of_flies_in_soc_space.csv           = "number_of_flies_in_soc_space",
  distance_traveled_between_interactions.csv = "distance_traveled_between_interactions",
  unique_partners_met_interaction_space.csv  = "unique_partners_met_interaction_space"
)

# Covariates that need a transform at read time (heavy right-skew): applied
# inside read_and_clean via this lookup. Only distance is log1p-transformed.
cov_transform <- c(
  distance_traveled_between_interactions.csv = "log1p"
)

# -----------------------------------------------------------------------------
# Model formula (fixed effects for the choice sub-model)
# Note: effect *names* in the formula are unchanged in the refactored goldfish.
# -----------------------------------------------------------------------------
# Decay analysis uses UNWEIGHTED inertia/recip only (see methods rationale):
# our comparison targets persistence (decay) of partner preference across
# treatments that may differ in baseline interaction rate; a count-weighted
# statistic would confound decay of social memory with between-treatment
# differences in interaction intensity. Unweighted isolates persistence of
# contact itself, is robust to heavy-tailed short-window counts, and gives a
# directly comparable memory-decay profile. Weighted is available as a
# robustness refit (choice_fixed_effects_weighted below), NOT in the main model.
# The Stage 0 diagnostic reports short-window repeat-count distributions so this
# choice can be confirmed empirically on the first run.
choice_fixed_effects <- dependent ~
  inertia(interaction_network, weighted = FALSE, window = 96) +
  inertia(interaction_network, weighted = FALSE, window = 288) +
  inertia(interaction_network, weighted = FALSE, window = 864) +
  inertia(interaction_network, weighted = FALSE, window = 2592) +

  recip(interaction_network, weighted = FALSE, window = 96) +
  recip(interaction_network, weighted = FALSE, window = 288) +
  recip(interaction_network, weighted = FALSE, window = 864) +
  recip(interaction_network, weighted = FALSE, window = 2592) +

  indeg(interaction_network, weighted = FALSE, window = 96) +
  indeg(interaction_network, weighted = FALSE, window = 288) +
  indeg(interaction_network, weighted = FALSE, window = 864) +
  indeg(interaction_network, weighted = FALSE, window = 2592) +

  outdeg(interaction_network, weighted = FALSE, window = 96) +
  outdeg(interaction_network, weighted = FALSE, window = 288) +
  outdeg(interaction_network, weighted = FALSE, window = 864) +
  outdeg(interaction_network, weighted = FALSE, window = 2592) +

  trans(interaction_network, window = 96) +
  trans(interaction_network, window = 288) +
  trans(interaction_network, window = 864) +
  trans(interaction_network, window = 2592) +

  # Covariate linking is now fixed (each covariate gets a distinct variable name
  # in the linking loop), so the five chosen exogenous covariates enter as alter
  # effects: do flies prefer partners high on these? distance is log1p at read.
  alter(nodesAttr$positive_influence) +
  alter(nodesAttr$negative_influence) +
  alter(nodesAttr$number_of_flies_in_soc_space) +
  alter(nodesAttr$distance_traveled_between_interactions) +
  alter(nodesAttr$unique_partners_met_interaction_space) +

  sim(nodesAttr$popularity)

# Robustness refit: identical model with WEIGHTED inertia/recip. Fit SEPARATELY
# (its own gather_model_data / estimate call), then compare the decay pattern
# against the unweighted fit. Do NOT combine both weightings in one model.
choice_fixed_effects_weighted <- dependent ~
  inertia(interaction_network, weighted = TRUE, window = 96) +
  inertia(interaction_network, weighted = TRUE, window = 288) +
  inertia(interaction_network, weighted = TRUE, window = 864) +
  inertia(interaction_network, weighted = TRUE, window = 2592) +

  recip(interaction_network, weighted = TRUE, window = 96) +
  recip(interaction_network, weighted = TRUE, window = 288) +
  recip(interaction_network, weighted = TRUE, window = 864) +
  recip(interaction_network, weighted = TRUE, window = 2592) +

  indeg(interaction_network, weighted = FALSE, window = 96) +
  indeg(interaction_network, weighted = FALSE, window = 288) +
  indeg(interaction_network, weighted = FALSE, window = 864) +
  indeg(interaction_network, weighted = FALSE, window = 2592) +

  outdeg(interaction_network, weighted = FALSE, window = 96) +
  outdeg(interaction_network, weighted = FALSE, window = 288) +
  outdeg(interaction_network, weighted = FALSE, window = 864) +
  outdeg(interaction_network, weighted = FALSE, window = 2592) +

  trans(interaction_network, window = 96) +
  trans(interaction_network, window = 288) +
  trans(interaction_network, window = 864) +
  trans(interaction_network, window = 2592) +

  # same five covariates as the unweighted formula
  alter(nodesAttr$positive_influence) +
  alter(nodesAttr$negative_influence) +
  alter(nodesAttr$number_of_flies_in_soc_space) +
  alter(nodesAttr$distance_traveled_between_interactions) +
  alter(nodesAttr$unique_partners_met_interaction_space) +

  sim(nodesAttr$popularity)

# -----------------------------------------------------------------------------
# RATE sub-model formula (ego-level).
# The rate model asks how fast a fly BECOMES ACTIVE, so effects are about the
# ego's own propensity — no alter/sim/recip/trans (those need a chosen partner).
# The leading `1 +` is the time intercept required by sub_model = "rate".
# Windows match the choice model (96/288/864/2592) so the two are comparable.
# Covariates that were constant/degenerate in the choice design are omitted
# here too; the rate rank-check (diagnose_rank on X_rate) will confirm before
# the full run. Add ego covariates back only after confirming they vary.
# -----------------------------------------------------------------------------
rate_fixed_effects <- dependent ~ 1 +
  indeg(interaction_network, weighted = FALSE, window = 96) +
  indeg(interaction_network, weighted = FALSE, window = 288) +
  indeg(interaction_network, weighted = FALSE, window = 864) +
  indeg(interaction_network, weighted = FALSE, window = 2592) +

  outdeg(interaction_network, weighted = FALSE, window = 96) +
  outdeg(interaction_network, weighted = FALSE, window = 288) +
  outdeg(interaction_network, weighted = FALSE, window = 864) +
  outdeg(interaction_network, weighted = FALSE, window = 2592) +

  # Same five exogenous covariates as choice, but ego-level here: does the fly's
  # OWN influence / spatial density / mobility / partner-diversity change how
  # fast it becomes active? distance is log1p at read time.
  ego(nodesAttr$positive_influence) +
  ego(nodesAttr$negative_influence) +
  ego(nodesAttr$number_of_flies_in_soc_space) +
  ego(nodesAttr$distance_traveled_between_interactions) +
  ego(nodesAttr$unique_partners_met_interaction_space)

# =============================================================================
# Helpers
# =============================================================================

# Build the actor/flies node data frame from an edgelist.
# Replaces the removed goldfish::getActors(): actor labels are the unique
# sender/receiver values, and all model attributes start at 0.
build_flies <- function(interaction_data, treatment_name,
                        treatment_flag = c(CsCh = "young", CS_10D = "old",
                                           Cs_5DIZ = "isolated")) {
  labels <- sort(unique(c(interaction_data$sender, interaction_data$receiver)))
  flies <- data.frame(label = labels, stringsAsFactors = FALSE)

  zero_attrs <- c(
    "popularity", "activity", "popularity_weighted", "activity_weighted",
    "positive_influence", "negative_influence",
    "positive_inf_weighted", "negative_inf_weighted",
    "distance_traveled_between_interactions", "number_of_flies_in_soc_space",
    "unique_partners_met_interaction_space", "unique_partners_met_social_space",
    "young", "old", "isolated"
  )
  flies[zero_attrs] <- 0L

  flag <- treatment_flag[[treatment_name]]
  if (!is.null(flag)) flies[[flag]] <- 1L
  flies
}

# Discover all edgelist CSVs across treatments.
list_group_files <- function(treatments, edgelist_path) {
  if (!dir.exists(edgelist_path)) {
    stop(glue::glue(
      "edgelist path does not exist:\n  {edgelist_path}\n",
      "Check cfg$base_path (currently points under it)."
    ), call. = FALSE)
  }

  rows <- lapply(treatments, function(tr) {
    tr_dir <- file.path(edgelist_path, tr)
    if (!dir.exists(tr_dir)) {
      warning(glue::glue("treatment folder not found: {tr_dir}"), call. = FALSE)
      return(NULL)
    }
    csvs <- list.files(tr_dir, pattern = "\\.csv$", full.names = TRUE)
    if (!length(csvs)) {
      warning(glue::glue("no CSVs in: {tr_dir}"), call. = FALSE)
      return(NULL)
    }
    data.frame(
      treatment_name = tr,
      file_path      = csvs,
      group_name     = tools::file_path_sans_ext(basename(csvs)),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    stop(glue::glue(
      "No edgelist CSVs found under {edgelist_path}.\n",
      "Looked in subfolders: {paste(treatments, collapse = ', ')}.\n",
      "Found there instead: {paste(list.files(edgelist_path), collapse = ', ')}\n",
      "-> Fix cfg$treatments to match the folder names, or point ",
      "cfg$base_path at the right directory."
    ), call. = FALSE)
  }
  out
}

# Preprocess a single group: build goldfish objects, link covariates,
# and run make_data() + gather_model_data() to produce the choice design matrix.
process_group <- function(i, files, cfg, paths, cov_attr_map,
                          choice_fixed_effects, log_file) {
  pid   <- Sys.getpid()
  t0    <- Sys.time()
  treatment_name <- files$treatment_name[i]
  group_name     <- files$group_name[i]

  cat(glue("{format(t0, '%Y-%m-%d %H:%M:%S')} [PID {pid}] START {group_name}"),
      "\n", file = log_file, append = TRUE, sep = "")

  # --- interaction events ---
  interaction_data <- read.csv(files$file_path[i])
  interaction_data <- interaction_data[, c("time", "sender", "receiver", "increment")]

  # --- nodes ---
  flies     <- build_flies(interaction_data, treatment_name)
  nodesAttr <- make_nodes(flies)

  # --- interaction network + dependent events ---
  interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
  interaction_network <- link_events(
    x = interaction_network, change_events = interaction_data, nodes = nodesAttr
  )
  dependent <- make_dependent_events(
    events = interaction_data, nodes = nodesAttr,
    default_network = interaction_network
  )

  # --- link node covariates ---
  group_cov_path <- file.path(paths$covariances, treatment_name, group_name)
  if (!dir.exists(group_cov_path)) {
    stop(glue(
      "covariance folder missing for group '{group_name}': {group_cov_path}\n",
      "Parent contains: ",
      "{paste(list.files(file.path(paths$covariances, treatment_name)), collapse = ', ')}"
    ), call. = FALSE)
  }
  for (fname in names(cov_attr_map)) {
    fpath <- file.path(group_cov_path, fname)
    if (!file.exists(fpath)) {
      stop(glue(
        "covariance file not found: {fpath}\n",
        "Folder actually contains: ",
        "{paste(list.files(group_cov_path), collapse = ', ')}"
      ), call. = FALSE)
    }
    # CRITICAL: goldfish's link_events tracks node change-events by the passed
    # VARIABLE NAME (substitute()). Reusing one name makes covariates 2..n look
    # like duplicates of the first and they are NOT applied (attrs stay at 0).
    # Assign each to a distinct name so every covariate links as its own events.
    attr_name <- cov_attr_map[[fname]]
    vname <- make.names(paste0("cov_", attr_name))
    assign(vname, read_and_clean(fname, group_cov_path))
    nodesAttr <- eval(bquote(
      link_events(nodesAttr, .(as.name(vname)), .(attr_name))))
  }

  # --- collect goldfish objects (native make_data collector) ---
  dataDynam <- make_data(dependent, interaction_network, nodesAttr, flies,
                         interaction_data)

  # --- preprocess the choice sub-model (native, fixed-effects only) ---
  # gather_model_data returns a tabular design: $stat_all_events (the design
  # matrix, columns named by $namesEffects), $selected, $n_candidates.
  group_data <- gather_model_data(
    formula   = choice_fixed_effects,
    model     = "DyNAM",
    sub_model = "choice",
    data      = dataDynam
  )

  t1 <- Sys.time()
  dur <- round(as.numeric(difftime(t1, t0, units = "secs")), 2)
  cat(glue("{format(t1, '%Y-%m-%d %H:%M:%S')} [PID {pid}] END   {group_name} ({dur}s)"),
      "\n", file = log_file, append = TRUE, sep = "")

  list(
    group_name     = group_name,
    treatment_name = treatment_name,
    group_index    = i,
    data_object    = group_data
  )
}

# -----------------------------------------------------------------------------
# make_cluster_safe(): create a PSOCK cluster, backing off if the OS refuses to
# fork that many processes ("Resource temporarily unavailable"). Tries the
# requested count, then progressively fewer, so a busy shared node degrades to a
# smaller cluster instead of erroring out.
# -----------------------------------------------------------------------------
make_cluster_safe <- function(n_cores) {
  tries <- unique(pmax(1L, c(n_cores, 8L, 4L, 2L, 1L)))
  tries <- tries[tries <= n_cores]
  for (n in tries) {
    cl <- tryCatch(parallel::makeCluster(n), error = function(e) NULL)
    if (!is.null(cl)) {
      if (n < n_cores) {
        message(glue("Requested {n_cores} cores but only {n} could be ",
                     "allocated (node busy / process limit). Using {n}."))
      }
      return(cl)
    }
    Sys.sleep(1)
  }
  stop("Could not create even a 1-worker cluster — the node is out of ",
       "process slots. Check `ulimit -u` and kill stray R workers ",
       "(pkill -u $USER -f parallel), then retry.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# diagnose_rank_design(): before a full refit, check ONE arena's freshly-built
# design (with the current formula + covariates) for constant or linearly
# dependent columns — the things that make estimation singular. Runs on one
# arena in seconds, so you catch collinearity (e.g. positive vs negative
# influence) BEFORE launching the 60-arena run.
#
# Rebuilds the arena design inline (same steps as process_group) so it reflects
# the CURRENT cov_attr_map and formula, not any stale checkpoint.
# -----------------------------------------------------------------------------
diagnose_rank_design <- function(cfg, paths, formula = choice_fixed_effects,
                                 sub_model = "choice", arena = NULL,
                                 cov_map = cov_attr_map) {
  files <- list_group_files(cfg$treatments, paths$edgelists)
  row <- if (is.null(arena)) files[1, ] else files[files$group_name == arena, ][1, ]
  g <- row$group_name; tr <- row$treatment_name
  ev <- read.csv(row$file_path)
  ev <- ev[, c("time", "sender", "receiver", "increment")]

  flies <- build_flies(ev, tr); nodesAttr <- make_nodes(flies)
  net <- make_network(nodes = nodesAttr, directed = TRUE)
  net <- link_events(x = net, change_events = ev, nodes = nodesAttr)
  dependent <- make_dependent_events(events = ev, nodes = nodesAttr,
                                     default_network = net)
  cov_dir <- file.path(paths$covariances, tr, g)
  for (fname in names(cov_map)) {
    attr_name <- cov_map[[fname]]; vname <- make.names(paste0("cov_", attr_name))
    assign(vname, read_and_clean(fname, cov_dir))
    nodesAttr <- eval(bquote(link_events(nodesAttr, .(as.name(vname)), .(attr_name))))
  }
  dd <- make_data(dependent, net, nodesAttr, flies, ev)
  X <- gather_model_data(formula, model = "DyNAM", sub_model = sub_model, data = dd)$stat_all_events

  v <- apply(X, 2, var)
  const <- names(v)[v == 0 | is.na(v)]
  Xnc <- X[, setdiff(colnames(X), const), drop = FALSE]
  Xsc <- scale(Xnc); Xsc[is.na(Xsc)] <- 0
  qrX <- qr(Xsc); dependent_cols <- colnames(Xnc)[qrX$pivot[-seq_len(qrX$rank)]]
  cm <- cor(Xnc); diag(cm) <- 0
  hi <- which(abs(cm) > 0.95, arr.ind = TRUE); hi <- hi[hi[,1] < hi[,2], , drop = FALSE]
  pairs <- if (nrow(hi)) data.frame(a = colnames(Xnc)[hi[,1]], b = colnames(Xnc)[hi[,2]],
    r = round(cm[hi], 3))[order(-abs(cm[hi])), ] else data.frame()

  cat(glue("\nArena {g} [{sub_model}]: {ncol(X)} cols, rank {qrX$rank} ",
           "=> deficiency {ncol(Xnc) - qrX$rank}\n\n"))
  cat("CONSTANT columns:\n"); print(const)
  cat("\nLINEARLY DEPENDENT (QR-flagged):\n"); print(dependent_cols)
  cat("\n|r| > 0.95 pairs:\n"); print(head(pairs, 20))
  invisible(list(constant = const, dependent = dependent_cols, pairs = pairs))
}

# =============================================================================
# STAGE 0 — Short-window repeat-count diagnostic (weighted vs unweighted)
# =============================================================================
# Decides empirically whether weighted and unweighted inertia/recip carry
# distinct information. For each window it reports, per directed dyad that has
# ANY interaction in that trailing window, the distribution of repeat counts:
#   - If counts are overwhelmingly 1 (esp. at window 96/288), weighted ≈
#     unweighted column-for-column -> keep UNWEIGHTED only (our default).
#   - If short-window counts are routinely >= 3, weighted carries independent
#     signal and keeping both in one model becomes defensible.
#
# Runs directly on the raw edgelists, no goldfish objects needed, so you can run
# it before Stage 1. Windows are in FRAMES (96=4s, 288=12s, 864=36s, 2592=108s).
run_repeat_count_diagnostic <- function(cfg, paths,
                                        windows = c(96, 288, 864, 2592)) {
  files <- list_group_files(cfg$treatments, paths$edgelists)
  if (!is.null(cfg$limit_groups)) files <- head(files, cfg$limit_groups)

  # For one group + one window, count events per (sender,receiver) that fall in
  # the trailing `w` frames of each event's time. We summarize the per-dyad
  # trailing counts pooled over all events (this is what a windowed weighted
  # statistic actually sees).
  per_group_window <- function(df, w) {
    ord <- order(df$time)
    tt  <- df$time[ord]; ss <- df$sender[ord]; rr <- df$receiver[ord]
    counts <- integer(0)
    for (k in seq_along(tt)) {
      in_win <- tt >= (tt[k] - w) & tt < tt[k] &
                ss == ss[k] & rr == rr[k]
      counts <- c(counts, sum(in_win))
    }
    counts <- counts[counts > 0]  # only dyads with any prior contact in-window
    if (!length(counts)) return(NULL)
    data.frame(
      window       = w,
      n_active     = length(counts),
      frac_exactly1 = mean(counts == 1),
      frac_ge3     = mean(counts >= 3),
      median_count = median(counts),
      p90_count    = as.numeric(quantile(counts, 0.90)),
      max_count    = max(counts)
    )
  }

  rows <- list()
  for (i in seq_len(nrow(files))) {
    df <- read.csv(files$file_path[i])
    df <- df[, c("time", "sender", "receiver")]
    for (w in windows) {
      s <- per_group_window(df, w)
      if (!is.null(s)) {
        s$group     <- files$group_name[i]
        s$treatment <- files$treatment_name[i]
        rows[[paste(files$group_name[i], w)]] <- s
      }
    }
  }
  detail <- do.call(rbind, rows)

  # Pooled summary per window across all groups (the decision table).
  agg <- do.call(rbind, lapply(split(detail, detail$window), function(d) {
    data.frame(
      window        = d$window[1],
      groups        = nrow(d),
      mean_frac_exactly1 = round(weighted.mean(d$frac_exactly1, d$n_active), 3),
      mean_frac_ge3      = round(weighted.mean(d$frac_ge3,      d$n_active), 3),
      median_of_medians  = median(d$median_count)
    )
  }))

  write.csv(detail, file.path(paths$tmp, "repeat_count_detail.csv"), row.names = FALSE)
  write.csv(agg,    file.path(paths$tmp, "repeat_count_summary.csv"), row.names = FALSE)
  message("Stage 0 complete. Decision table:")
  print(agg)
  message(
    "Rule of thumb: if mean_frac_exactly1 is high (>~0.8) at windows 96/288, ",
    "weighted ~ unweighted -> keep UNWEIGHTED only. If mean_frac_ge3 is ",
    "non-trivial (>~0.2), weighted carries independent signal."
  )
  invisible(list(detail = detail, summary = agg))
}

# =============================================================================
# STAGE 1 — Parallel per-group preprocessing
# =============================================================================
run_preprocessing <- function(cfg, paths,
                              formula  = choice_fixed_effects,
                              log_file = file.path(paths$results,
                                                   "parallel_execution_log.txt"),
                              cov_map  = cov_attr_map) {
  files <- list_group_files(cfg$treatments, paths$edgelists)
  if (!is.null(cfg$limit_groups)) files <- head(files, cfg$limit_groups)
  message(glue("Preprocessing {nrow(files)} groups on {cfg$n_cores} cores."))
  message(glue("Log: {log_file}"))

  cl <- make_cluster_safe(cfg$n_cores)
  on.exit({ stopCluster(cl); registerDoSEQ() }, add = TRUE)
  registerDoParallel(cl)

  parallel_results <- foreach(
    i = seq_len(nrow(files)),
    .packages     = c("goldfish", "glue", "tools"),
    .export        = c("process_group", "build_flies", "read_and_clean"),
    .errorhandling = "pass"
  ) %dopar% {
    # formula, log_file, cov_map are LOCALS of run_preprocessing, captured
    # automatically by foreach. build_flies carries its treatment map as a
    # default arg, so no global lookup is needed on the worker.
    process_group(i, files, cfg, paths, cov_map, formula, log_file)
  }

  # Split successes from errors.
  data_list <- list()
  info_rows <- list()
  n_failed  <- 0L
  first_err <- NULL
  for (res in parallel_results) {
    if (inherits(res, "list") && !is.null(res$group_name)) {
      g <- res$group_name
      data_list[[g]] <- res$data_object
      info_rows[[g]] <- data.frame(
        group   = g,
        treatment = res$treatment_name,
        ixGroup = res$group_index,
        stringsAsFactors = FALSE
      )
    } else {
      n_failed <- n_failed + 1L
      if (is.null(first_err)) first_err <- conditionMessage(res)
    }
  }
  if (n_failed > 0) {
    message(glue("{n_failed} group(s) failed. First error (representative):"))
    message(first_err)
  }
  info <- do.call(rbind, info_rows)

  n_ok <- length(data_list)
  if (n_ok == 0) {
    stop(glue(
      "Stage 1: ALL {nrow(files)} groups failed — nothing to save. ",
      "Scroll up for the first 'Group failed:' message (now includes the ",
      "missing file/folder). Fix that before rerunning downstream stages."
    ), call. = FALSE)
  }
  if (n_ok < nrow(files)) {
    message(glue("Stage 1: {n_ok}/{nrow(files)} groups succeeded ",
                 "({nrow(files) - n_ok} failed)."))
  }

  saveRDS(data_list, file.path(paths$tmp, "all_groups_data_list.rds"))
  saveRDS(info,      file.path(paths$tmp, "all_groups_info.rds"))
  message("Stage 1 complete. Saved to ", paths$tmp)

  list(data_list = data_list, info = info)
}

# =============================================================================
# STAGE 2 — Collinearity diagnostics across groups
# =============================================================================
run_collinearity_plot <- function(paths, cor_thresh) {
  data_list <- readRDS(file.path(paths$tmp, "all_groups_data_list.rds"))

  # Native gather_model_data output: design matrix is $stat_all_events, with
  # columns already named by the effects. Choice sub-model has no intercept
  # column, so nothing to strip.
  design_of <- function(g) g$stat_all_events

  all_cols <- sort(unique(unlist(lapply(data_list, \(g) colnames(design_of(g))))))
  count_mat <- matrix(0, length(all_cols), length(all_cols),
                      dimnames = list(all_cols, all_cols))

  for (g in data_list) {
    X <- design_of(g)
    const <- which(apply(X, 2, var) == 0)
    if (length(const)) X <- X[, -const, drop = FALSE]
    if (ncol(X) < 2) next
    adj <- (abs(cor(X)) > cor_thresh) * 1
    diag(adj) <- 0
    p <- colnames(adj)
    count_mat[p, p] <- count_mat[p, p] + adj
  }

  keep <- (rowSums(count_mat) > 0) | (colSums(count_mat) > 0)
  to_plot <- count_mat[keep, keep, drop = FALSE]

  out_png <- file.path(paths$tmp, "high_correlation_plot.png")
  png(out_png, width = 14, height = 14, units = "in", res = 300)
  corrplot::corrplot(
    to_plot, method = "color", type = "upper",
    col = colorRampPalette(c("white", "red"))(100),
    is.corr = FALSE, addgrid.col = "grey",
    tl.col = "black", tl.srt = 90, tl.cex = 0.9,
    title = as.character(cor_thresh), mar = c(0, 0, 1, 0)
  )
  dev.off()
  message("Stage 2 complete: ", out_png)
  invisible(count_mat)
}

# =============================================================================
# STAGE 3 — Native DyNAM choice estimation (fixed effects)
# =============================================================================
# On this branch the estimator is goldfish::estimate_dynam(); no Cox / Stan
# layer. We estimate per group and then (optionally) a pooled fit. Because the
# formula references objects (interaction_network, nodesAttr, flies) that only
# exist inside process_group(), estimation is cleanest done in the SAME env the
# objects are built. So Stage 3 re-runs the object construction per group and
# calls estimate_dynam() directly, rather than reloading the gathered matrices.
#
# `formula` lets you pass choice_fixed_effects (main, unweighted) or
# choice_fixed_effects_weighted (robustness refit). Run it twice and compare the
# decay pattern across the two fits — this is the valid way to "test both"
# weightings without putting collinear terms in one model.
run_estimation <- function(cfg, paths, formula, tag = "unweighted",
                           sub_model = "choice",
                           cov_map = cov_attr_map) {
  files <- list_group_files(cfg$treatments, paths$edgelists)
  if (!is.null(cfg$limit_groups)) files <- head(files, cfg$limit_groups)

  cl <- make_cluster_safe(cfg$n_cores)
  on.exit({ stopCluster(cl); registerDoSEQ() }, add = TRUE)
  registerDoParallel(cl)

  fits <- foreach(
    i = seq_len(nrow(files)),
    .packages = c("goldfish", "glue", "tools"),
    .export        = c("build_flies", "read_and_clean"),
    .errorhandling = "pass"
  ) %dopar% {

    treatment_name <- files$treatment_name[i]
    group_name     <- files$group_name[i]

    interaction_data <- read.csv(files$file_path[i])
    interaction_data <- interaction_data[, c("time", "sender", "receiver", "increment")]

    flies     <- build_flies(interaction_data, treatment_name)
    nodesAttr <- make_nodes(flies)
    interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
    interaction_network <- link_events(
      x = interaction_network, change_events = interaction_data, nodes = nodesAttr)
    dependent <- make_dependent_events(
      events = interaction_data, nodes = nodesAttr,
      default_network = interaction_network)

    group_cov_path <- file.path(paths$covariances, treatment_name, group_name)
    for (fname in names(cov_map)) {
      # distinct variable name per covariate — see process_group for why.
      attr_name <- cov_map[[fname]]
      vname <- make.names(paste0("cov_", attr_name))
      assign(vname, read_and_clean(fname, group_cov_path))
      nodesAttr <- eval(bquote(
        link_events(nodesAttr, .(as.name(vname)), .(attr_name))))
    }

    dataDynam <- make_data(dependent, interaction_network, nodesAttr, flies,
                           interaction_data)

    est <- estimate_dynam(
      x         = formula,
      sub_model = sub_model,
      data      = dataDynam
    )
    list(group_name = group_name, treatment_name = treatment_name, fit = est)
  }

  # Collect coefficient tables (broom::tidy works on goldfish result objects).
  coef_rows <- list()
  n_fail <- 0L; first_err <- NULL
  for (res in fits) {
    if (inherits(res, "list") && !is.null(res$fit)) {
      tb <- tryCatch(broom::tidy(res$fit), error = function(e) NULL)
      if (!is.null(tb)) {
        tb$group     <- res$group_name
        tb$treatment <- res$treatment_name
        coef_rows[[res$group_name]] <- tb
      } else {
        n_fail <- n_fail + 1L
        if (is.null(first_err)) first_err <- "broom::tidy() returned NULL for a fit"
      }
    } else {
      n_fail <- n_fail + 1L
      if (is.null(first_err)) first_err <- conditionMessage(res)
    }
  }
  if (length(coef_rows) == 0) {
    stop(glue(
      "Stage 3 [{sub_model}]: ALL groups failed — no coefficients to write.\n",
      "First error (representative): {first_err}"
    ), call. = FALSE)
  }
  if (n_fail > 0) {
    message(glue("Stage 3 [{sub_model}]: {length(coef_rows)} succeeded, ",
                 "{n_fail} failed. First error: {first_err}"))
  }
  coefs <- do.call(rbind, coef_rows)
  # File naming: choice keeps the weighting tag; rate uses sub_model in the name.
  stem <- if (sub_model == "choice") glue("choice_coefs_{tag}") else glue("rate_coefs_{tag}")
  out_csv <- file.path(paths$tmp, glue("{stem}.csv"))
  write.csv(coefs, out_csv, row.names = FALSE)
  saveRDS(fits, file.path(paths$tmp, glue("{stem}_fits.rds")))
  message(glue("Stage 3 [{sub_model}] ({tag}) complete: {out_csv}"))
  invisible(coefs)
}

# =============================================================================
# STAGE 4 — Simple meta-analysis across arenas (per treatment)
# =============================================================================
# Per mentor's guidance: fit one model per ARENA (done in Stage 3 — each group
# is an arena), then combine the per-arena coefficients with a simple,
# closed-form meta-analysis. NO Stan / brms — this is classic inverse-variance
# pooling (DerSimonian-Laird random effects), runs in milliseconds.
#
# For each term, WITHIN each treatment, we pool the arena-level estimates:
#   - fixed-effect weight  w_i   = 1 / SE_i^2
#   - Q (heterogeneity), tau^2 (DL between-arena variance)
#   - random-effects weight w*_i = 1 / (SE_i^2 + tau^2)
#   - pooled estimate, SE, 95% CI, and I^2 (share of variance from heterogeneity)
#
# Reads the Stage 3 CSV (choice_coefs_<tag>.csv), so it is instant and needs no
# re-estimation.

# DerSimonian-Laird random-effects pool for one set of (estimate, se) values.
meta_pool_dl <- function(est, se) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]; se <- se[ok]
  k <- length(est)
  if (k == 0) return(NULL)
  if (k == 1) {
    return(data.frame(k = 1, pooled = est, se = se,
                      ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
                      Q = NA_real_, tau2 = 0, I2 = NA_real_))
  }
  w   <- 1 / se^2
  fe  <- sum(w * est) / sum(w)                       # fixed-effect mean
  Q   <- sum(w * (est - fe)^2)                       # Cochran's Q
  df  <- k - 1
  C   <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)                       # DL between-study variance
  wr  <- 1 / (se^2 + tau2)                           # random-effects weights
  pooled <- sum(wr * est) / sum(wr)
  se_p   <- sqrt(1 / sum(wr))
  I2 <- if (Q > df) (Q - df) / Q * 100 else 0        # % variance from heterogeneity
  data.frame(k = k, pooled = pooled, se = se_p,
             ci_low = pooled - 1.96 * se_p, ci_high = pooled + 1.96 * se_p,
             Q = Q, tau2 = tau2, I2 = I2)
}

# Parse "interaction_network indeg 2592" -> effect="indeg", window=2592.
# Non-windowed terms (e.g. "nodesAttr popularity sim") get window = NA.
parse_term <- function(term) {
  win <- suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", term)))
  eff <- term
  eff <- sub("^interaction_network\\s+", "", eff)
  eff <- sub("\\s+\\d+$", "", eff)                   # strip trailing window
  data.frame(effect = eff, window = win, stringsAsFactors = FALSE)
}

run_meta_analysis <- function(paths, tag = "unweighted", sub_model = "choice") {
  stem <- if (sub_model == "choice") glue("choice_coefs_{tag}") else glue("rate_coefs_{tag}")
  in_csv <- file.path(paths$tmp, glue("{stem}.csv"))
  if (!file.exists(in_csv)) {
    stop("No Stage 3 output at ", in_csv,
         " — run run_estimation(..., sub_model = '", sub_model,
         "', tag = '", tag, "') first.", call. = FALSE)
  }
  co <- read.csv(in_csv, stringsAsFactors = FALSE)

  # Pool within each (treatment, term).
  keys <- unique(co[, c("treatment", "term")])
  rows <- lapply(seq_len(nrow(keys)), function(j) {
    tr <- keys$treatment[j]; tm <- keys$term[j]
    sub <- co[co$treatment == tr & co$term == tm, ]
    p <- meta_pool_dl(sub$estimate, sub$std.error)
    if (is.null(p)) return(NULL)
    cbind(data.frame(treatment = tr, term = tm), parse_term(tm), p)
  })
  meta <- do.call(rbind, rows)
  meta <- meta[order(meta$effect, meta$window, meta$treatment), ]

  out_csv <- file.path(paths$tmp, glue("meta_pooled_{sub_model}_{tag}.csv"))
  write.csv(meta, out_csv, row.names = FALSE)
  message(glue("Stage 4 [{sub_model}] ({tag}) complete: {out_csv}"))

  # Decay figure: pooled coefficient vs window, one line per treatment, one
  # panel per windowed effect.
  fig <- tryCatch({
    library(ggplot2)
    dfw <- meta[is.finite(meta$window), ]
    p <- ggplot(dfw, aes(window, pooled, colour = treatment, group = treatment)) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = treatment),
                  alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
      scale_x_log10(breaks = c(96, 288, 864, 2592)) +
      facet_wrap(~ effect, scales = "free_y") +
      labs(x = "window (frames, log scale)", y = "pooled coefficient",
           title = glue("Decay of {sub_model} effects by treatment ({tag})"),
           subtitle = "arena-level DyNAM fits, DL random-effects pooled") +
      theme_bw(base_size = 11)
    out_png <- file.path(paths$tmp, glue("decay_curves_{sub_model}_{tag}.png"))
    ggsave(out_png, p, width = 11, height = 8, dpi = 200)
    message(glue("Stage 4 decay figure: {out_png}"))
    out_png
  }, error = function(e) {
    message("Decay plot skipped (", conditionMessage(e), ") — CSV still written.")
    NULL
  })

  # Covariate figure: NON-windowed effects (covariates, intercept) get a forest
  # plot instead — pooled estimate + 95% CI per effect, one row per treatment.
  # A decay curve makes no sense for these, but they still need visualizing.
  fig_cov <- tryCatch({
    library(ggplot2)
    dfc <- meta[!is.finite(meta$window), ]
    if (nrow(dfc) == 0) {
      message("No non-windowed (covariate) effects to plot.")
      NULL
    } else {
      p2 <- ggplot(dfc, aes(x = pooled, y = effect, colour = treatment)) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                       height = 0.25, position = position_dodge(width = 0.6)) +
        geom_point(size = 2, position = position_dodge(width = 0.6)) +
        labs(x = "pooled coefficient (95% CI)", y = NULL,
             title = glue("Covariate effects by treatment — {sub_model} ({tag})"),
             subtitle = "arena-level DyNAM fits, DL random-effects pooled") +
        theme_bw(base_size = 11)
      out_png2 <- file.path(paths$tmp, glue("covariate_effects_{sub_model}_{tag}.png"))
      ggsave(out_png2, p2, width = 9, height = 6, dpi = 200)
      message(glue("Stage 4 covariate figure: {out_png2}"))
      out_png2
    }
  }, error = function(e) {
    message("Covariate plot skipped (", conditionMessage(e), ") — CSV still written.")
    NULL
  })

  invisible(meta)
}

# =============================================================================
# Checkpointing
# =============================================================================
# Stage 1 writes two checkpoints to paths$tmp:
#   all_groups_data_list.rds  — per-group design matrices
#   all_groups_info.rds       — group/treatment index
# Stage 2 and (optionally) later analysis read these back, so Stage 1 only needs
# to run ONCE. Sourcing this file does NOT run anything — it just loads the
# functions, cfg, and paths. Call the stage functions yourself.

# Are the Stage 1 checkpoints already on disk?
stage1_done <- function(paths) {
  file.exists(file.path(paths$tmp, "all_groups_data_list.rds")) &&
    file.exists(file.path(paths$tmp, "all_groups_info.rds"))
}

# Load the Stage 1 checkpoints into a list (for interactive inspection).
load_stage1 <- function(paths) {
  if (!stage1_done(paths)) {
    stop("Stage 1 checkpoints not found in ", paths$tmp,
         " — run run_preprocessing(cfg, paths) first.", call. = FALSE)
  }
  list(
    data_list = readRDS(file.path(paths$tmp, "all_groups_data_list.rds")),
    info      = readRDS(file.path(paths$tmp, "all_groups_info.rds"))
  )
}

# =============================================================================
# Run helpers — call these explicitly; sourcing the file runs nothing.
# =============================================================================

# Full pipeline from scratch. Skips Stage 1 automatically if its checkpoints
# already exist (set force_preprocess = TRUE to recompute them anyway).
run_all <- function(cfg, paths, force_preprocess = FALSE) {
  run_repeat_count_diagnostic(cfg, paths)                 # Stage 0
  if (force_preprocess || !stage1_done(paths)) {
    run_preprocessing(cfg, paths)                         # Stage 1 (expensive)
  } else {
    message("Stage 1 checkpoint found — skipping preprocessing. ",
            "Pass force_preprocess = TRUE to recompute.")
  }
  run_collinearity_plot(paths, cfg$cor_thresh)            # Stage 2
  run_estimation(cfg, paths, choice_fixed_effects, tag = "unweighted")  # Stage 3
  run_meta_analysis(paths, tag = "unweighted")            # Stage 4
  invisible(TRUE)
}

# Continue from Stage 2 onward, reusing the Stage 1 checkpoint. This is the
# "next session" entry point:
#   source("<this file>"); resume_from_stage2(cfg, paths)
resume_from_stage2 <- function(cfg, paths) {
  if (!stage1_done(paths)) {
    stop("No Stage 1 checkpoint in ", paths$tmp,
         " — nothing to resume. Run run_preprocessing(cfg, paths) first.",
         call. = FALSE)
  }
  run_collinearity_plot(paths, cfg$cor_thresh)            # Stage 2
  run_estimation(cfg, paths, choice_fixed_effects, tag = "unweighted")  # Stage 3
  run_meta_analysis(paths, tag = "unweighted")            # Stage 4
  invisible(TRUE)
}

# NOTE: this file intentionally does NOT auto-run on source(). To execute:
#   run_all(cfg, paths)            # first time, from scratch (choice)
#   resume_from_stage2(cfg, paths) # later sessions, reuse Stage 1 checkpoint
#
# CHOICE model:
#   run_estimation(cfg, paths, choice_fixed_effects, tag = "unweighted")
#   run_estimation(cfg, paths, choice_fixed_effects_weighted, tag = "weighted")
#   run_meta_analysis(paths, tag = "unweighted")   # decay curves + pooled CSV
#   run_meta_analysis(paths, tag = "weighted")
#
# RATE model (Option A: per-arena estimate_dynam + same DL meta-analysis):
#   run_estimation(cfg, paths, rate_fixed_effects, tag = "main",
#                  sub_model = "rate")
#   run_meta_analysis(paths, tag = "main", sub_model = "rate")
#   # -> rate_coefs_main.csv, meta_pooled_rate_main.csv, decay_curves_rate_main.png