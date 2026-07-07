# =============================================================================
# STANDALONE covariate diagnostic
#
# Answers: your covariates vary in the raw files and are correctly time-aligned,
# yet some came through the choice design as CONSTANT (which forced us to drop
# them). This script settles whether they are *really* constant in the scored
# design, or whether the earlier "constant" flag came from a stale checkpoint.
#
# It is fully SELF-CONTAINED: no need to source the main pipeline. Just set the
# CONFIG block below and run. Requires only the goldfish package
# (@refactor/rate_prep) and glue.
#
# Run:
#   source("/home/milky/droso-pipe/6_dynam_models/src/cov_diagnostic_standalone.R")
#   # then either / both:
#   diag_raw   <- check_raw_and_time()      # raw variation + time overlap
#   diag_design <- check_design_variance()  # DECISIVE: constant in scored design?
# =============================================================================

suppressPackageStartupMessages({
  library(goldfish)
  library(glue)
})

# -----------------------------------------------------------------------------
# CONFIG — edit these to match your machine.
# -----------------------------------------------------------------------------
BASE_PATH   <- "/home/milky/droso-pipe/6_dynam_models"
TREATMENTS  <- c("CS_10D", "Cs_5DIZ", "CsCh")
EDGELIST_DIR <- file.path(BASE_PATH, "data", "edgelists")
COV_DIR_ROOT <- file.path(BASE_PATH, "data", "covariances")

# Treatment -> experimental flag (for build_flies).
TREATMENT_FLAG <- c(CsCh = "young", CS_10D = "old", Cs_5DIZ = "isolated")

# Full covariate map: file -> attribute name (all 11).
COV_MAP <- c(
  in_degree.csv                              = "popularity",
  out_degree.csv                             = "activity",
  in_weighted.csv                            = "popularity_weighted",
  out_weighted.csv                           = "activity_weighted",
  positiveinfluence.csv                      = "positive_influence",
  negativeinfluence.csv                      = "negative_influence",
  positiveinf_weighted.csv                   = "positive_inf_weighted",
  negativeinf_weighted.csv                   = "negative_inf_weighted",
  number_of_flies_in_soc_space.csv           = "number_of_flies_in_soc_space",
  unique_partners_met_interaction_space.csv  = "unique_partners_met_interaction_space",
  unique_partners_met_social_space.csv       = "unique_partners_met_social_space"
)

# -----------------------------------------------------------------------------
# Self-contained helpers (copied from the pipeline so this file stands alone).
# -----------------------------------------------------------------------------

# Load a covariate change-event CSV (time, actor/node, value) into goldfish's
# expected shape: time, node, replace.
.read_and_clean <- function(filename, directory) {
  df <- read.csv(file.path(directory, filename), stringsAsFactors = FALSE)
  if ("actor" %in% names(df)) names(df)[names(df) == "actor"] <- "node"
  if ("value" %in% names(df)) names(df)[names(df) == "value"] <- "replace"
  df <- df[, c("time", "node", "replace")]
  df$time    <- as.numeric(df$time)
  df$replace <- as.numeric(df$replace)
  df[order(df$time), ]
}

# Build the flies node data frame with all covariate attributes at 0.
.build_flies <- function(interaction_data, treatment_name) {
  labels <- sort(unique(c(interaction_data$sender, interaction_data$receiver)))
  flies <- data.frame(label = labels, stringsAsFactors = FALSE)
  flies[unname(COV_MAP)] <- 0
  flies[c("young", "old", "isolated")] <- 0L
  flag <- TREATMENT_FLAG[[treatment_name]]
  if (!is.null(flag)) flies[[flag]] <- 1L
  flies
}

# Discover edgelist CSVs across treatments.
.list_group_files <- function() {
  rows <- lapply(TREATMENTS, function(tr) {
    csvs <- list.files(file.path(EDGELIST_DIR, tr),
                       pattern = "\\.csv$", full.names = TRUE)
    if (!length(csvs)) return(NULL)
    data.frame(treatment_name = tr, file_path = csvs,
               group_name = tools::file_path_sans_ext(basename(csvs)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Recompute the scored design for one event sequence via goldfish's own
# gather_model_data (the trusted, estimation-identical path).
.recompute_design <- function(events, formula, treatment_name, cov_dir) {
  events <- events[, c("time", "sender", "receiver", "increment")]
  flies     <- .build_flies(events, treatment_name)
  nodesAttr <- make_nodes(flies)
  interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
  interaction_network <- link_events(
    x = interaction_network, change_events = events, nodes = nodesAttr)
  dependent <- make_dependent_events(
    events = events, nodes = nodesAttr, default_network = interaction_network)

  # CRITICAL FIX: goldfish's link_events tracks node change-events by the NAME
  # of the variable passed (via substitute()). Reusing one name ("cov_df") for
  # all covariates makes goldfish treat covariates 2..n as duplicates of the
  # first and silently NOT apply them -> those attributes stay at init 0
  # (constant). Assign each covariate to a DISTINCT variable name so each is
  # linked as its own event set.
  for (fname in names(COV_MAP)) {
    attr_name <- COV_MAP[[fname]]
    vname <- make.names(paste0("cov_", attr_name))
    assign(vname, .read_and_clean(fname, cov_dir))
    call_link <- bquote(
      link_events(nodesAttr, .(as.name(vname)), .(attr_name)))
    nodesAttr <- eval(call_link)
  }

  dataDynam <- make_data(dependent, interaction_network, nodesAttr, flies, events)
  gmd <- gather_model_data(formula = formula, model = "DyNAM",
                           sub_model = "choice", data = dataDynam)
  gmd$stat_all_events
}

# Pick the arena to inspect: NULL -> first arena found.
.pick_arena <- function(arena) {
  files <- .list_group_files()
  if (is.null(files)) stop("No edgelists found under ", EDGELIST_DIR, call. = FALSE)
  if (is.null(arena)) files[1, ] else {
    r <- files[files$group_name == arena, ]
    if (!nrow(r)) stop("Arena not found: ", arena, call. = FALSE)
    r[1, ]
  }
}

# =============================================================================
# CHECK 1 — raw variation + time overlap (no modeling, just file inspection).
# =============================================================================
check_raw_and_time <- function(arena = NULL) {
  row <- .pick_arena(arena)
  g <- row$group_name; tr <- row$treatment_name
  cov_dir <- file.path(COV_DIR_ROOT, tr, g)
  ev <- read.csv(row$file_path)
  ev_tmin <- min(ev$time); ev_tmax <- max(ev$time)

  message(glue("Arena: {g} ({tr}) | events n={nrow(ev)} time [{ev_tmin}, {ev_tmax}]"))

  out <- data.frame()
  for (fname in names(COV_MAP)) {
    fpath <- file.path(cov_dir, fname)
    if (!file.exists(fpath)) {
      out <- rbind(out, data.frame(file = fname, attr = COV_MAP[[fname]],
        raw_varies = NA, raw_range = "MISSING", cov_tmin = NA, cov_tmax = NA,
        time_overlaps = NA, stringsAsFactors = FALSE)); next
    }
    cdf <- read.csv(fpath)
    vcol <- if ("value" %in% names(cdf)) "value" else tail(names(cdf), 1)
    tcol <- if ("time"  %in% names(cdf)) "time"  else names(cdf)[1]
    vr <- range(cdf[[vcol]], na.rm = TRUE); trg <- range(cdf[[tcol]], na.rm = TRUE)
    out <- rbind(out, data.frame(file = fname, attr = COV_MAP[[fname]],
      raw_varies = diff(vr) > 0, raw_range = glue("[{round(vr[1],2)}, {round(vr[2],2)}]"),
      cov_tmin = trg[1], cov_tmax = trg[2],
      time_overlaps = !(trg[2] < ev_tmin | trg[1] > ev_tmax),
      stringsAsFactors = FALSE))
  }
  cat("\n===== RAW FILE + TIME-OVERLAP =====\n"); print(out, row.names = FALSE)
  invisible(out)
}

# =============================================================================
# CHECK 2 — DECISIVE: are the covariate effects constant in the SCORED design?
# =============================================================================
check_design_variance <- function(arena = NULL,
                                  effect_type = c("sim", "alter", "ego")) {
  effect_type <- match.arg(effect_type)
  row <- .pick_arena(arena)
  g <- row$group_name; tr <- row$treatment_name
  cov_dir <- file.path(COV_DIR_ROOT, tr, g)
  ev <- read.csv(row$file_path)

  terms <- vapply(unname(COV_MAP),
                  function(a) glue("{effect_type}(nodesAttr${a})"), character(1))
  fml <- as.formula(paste("dependent ~", paste(terms, collapse = " + ")))

  message(glue("Recomputing design for {g} with {length(terms)} '{effect_type}' ",
               "covariate effects (this scores through goldfish)..."))
  X <- .recompute_design(ev, fml, tr, cov_dir)

  v <- apply(X, 2, var)
  res <- data.frame(column = colnames(X), variance = round(v, 6),
    status = ifelse(v == 0 | is.na(v), "CONSTANT", "varies"),
    stringsAsFactors = FALSE)
  cat("\n===== DESIGN-LEVEL VARIANCE (scored via goldfish) =====\n")
  print(res, row.names = FALSE)
  n_const <- sum(res$status == "CONSTANT")
  cat(glue("\n{n_const}/{nrow(res)} covariate effects CONSTANT in the scored design.\n"))
  if (n_const == 0) {
    cat("=> None constant. The earlier drop was based on a stale checkpoint.\n",
        "   RESTORE these covariates to the model.\n", sep = "")
  } else {
    cat("=> These are genuinely flat in the design despite varying raw data:\n", sep = "")
    print(res[res$status == "CONSTANT", "column"])
  }
  invisible(res)
}

message("Standalone covariate diagnostic loaded. Run:")
message("  check_raw_and_time()        # raw variation + time overlap")
message("  check_design_variance()     # DECISIVE: constant in scored design?")
message("  check_design_variance(effect_type = 'alter')  # try alter/ego too")
message("  time_covariate_linking()    # times each covariate ADDED one-by-one")

# =============================================================================
# time_covariate_linking(): add covariates ONE AT A TIME and time each
# recompute, using a SMALL slice of events so it returns quickly. Distinguishes
# a true hang (one covariate never returns) from genuine per-covariate cost.
# max_events caps the arena to keep each step fast.
# =============================================================================
time_covariate_linking <- function(arena = NULL, effect_type = "sim",
                                    max_events = 300L, n_cov = length(COV_MAP)) {
  row <- .pick_arena(arena)
  g <- row$group_name; tr <- row$treatment_name
  cov_dir <- file.path(COV_DIR_ROOT, tr, g)
  ev <- read.csv(row$file_path)
  if (nrow(ev) > max_events) ev <- ev[seq_len(max_events), ]
  message(glue("Arena {g}: timing with {nrow(ev)} events, up to {n_cov} covariates."))

  cov_names <- names(COV_MAP)[seq_len(min(n_cov, length(COV_MAP)))]
  for (k in seq_along(cov_names)) {
    use <- cov_names[seq_len(k)]
    terms <- vapply(unname(COV_MAP[use]),
                    function(a) glue("{effect_type}(nodesAttr${a})"), character(1))
    fml <- as.formula(paste("dependent ~", paste(terms, collapse = " + ")))
    t0 <- Sys.time()
    res <- tryCatch({
      # inline recompute restricted to `use` covariates
      evs <- ev[, c("time", "sender", "receiver", "increment")]
      flies <- .build_flies(evs, tr); nodesAttr <- make_nodes(flies)
      net <- make_network(nodes = nodesAttr, directed = TRUE)
      net <- link_events(x = net, change_events = evs, nodes = nodesAttr)
      dep <- make_dependent_events(events = evs, nodes = nodesAttr,
                                   default_network = net)
      for (fn in use) {
        an <- COV_MAP[[fn]]; vn <- make.names(paste0("cov_", an))
        assign(vn, .read_and_clean(fn, cov_dir))
        nodesAttr <- eval(bquote(link_events(nodesAttr, .(as.name(vn)), .(an))))
      }
      dd <- make_data(dep, net, nodesAttr, flies, evs)
      X <- gather_model_data(fml, model = "DyNAM", sub_model = "choice", data = dd)$stat_all_events
      "OK"
    }, error = function(e) paste("ERR:", conditionMessage(e)))
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    message(glue("  {k} covariate(s) [+{tail(use,1)}]: {dt}s  {res}"))
  }
  invisible(NULL)
}