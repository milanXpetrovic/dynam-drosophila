# =============================================================================
# DyNAM: Fly interaction network — model fitting + GoF via forward simulation
# Data: CsCh_A1_24_05_2023-12_27.csv
# Models:
#   Rate:   intercept + indeg(window=200) + outdeg(window=200)
#   Choice: inertia(window=200) + indeg(window=200) + outdeg(window=200)
# GoF:    100 simulated networks vs. observed on:
#           indegree dist, outdegree dist, tie weight dist,
#           reciprocity dist, weighted density
# =============================================================================

library(goldfish)
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(2012)

data_path <- "C:/Complex Networks/flies/simulation attempt/CsCh_A1_24_05_2023-12_27.csv"
out_dir   <- "C:/Complex Networks/flies/simulation attempt"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

W <- 200L   # window size (in timestamp units)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

raw <- read.csv(data_path, stringsAsFactors = FALSE)

# Clean column names — drop row index column if present
if ("X" %in% names(raw)) raw <- raw[, !names(raw) %in% "X"]

# Ensure correct types
events_df <- raw %>%
  transmute(
    time      = as.numeric(time),
    sender    = as.character(sender),
    receiver  = as.character(receiver),
    increment = as.integer(increment)
  ) %>%
  arrange(time)

cat("Events loaded:", nrow(events_df), "\n")
cat("Time range: [", min(events_df$time), ",", max(events_df$time), "]\n")

# All unique actors — discovered from data, not hardcoded
all_actors <- sort(unique(c(events_df$sender, events_df$receiver)))
N          <- length(all_actors)
cat("Unique actors:", N, ":", paste(all_actors, collapse = ", "), "\n")

T_max <- max(events_df$time)


# =============================================================================
# 2. DEFINE GOLDFISH OBJECTS
# =============================================================================

actors_df <- data.frame(
  label   = all_actors,
  present = TRUE,
  stringsAsFactors = FALSE
)

actors_gf <- defineNodes(actors_df)

fly_network <- defineNetwork(nodes = actors_gf, directed = TRUE) %>%
  linkEvents(changeEvents = events_df, nodes = actors_gf)

fly_dependent <- defineDependentEvents(
  events         = events_df,
  nodes          = actors_gf,
  defaultNetwork = fly_network
)


# =============================================================================
# 3. FIT DYNAM MODELS
# =============================================================================

cat("\n--- Fitting rate model ---\n")
mod_rate <- estimate(
  fly_dependent ~ 1 +
    indeg(fly_network,  window = W) +
    outdeg(fly_network, window = W),
  model    = "DyNAM",
  subModel = "rate",
  estimationInit = list(engine = "default_c")
)
print(summary(mod_rate))

cat("\n--- Fitting choice model ---\n")
mod_choice <- estimate(
  fly_dependent ~
    inertia(fly_network, window = W) +
    indeg(fly_network,   window = W) +
    outdeg(fly_network,  window = W),
  model    = "DyNAM",
  subModel = "choice",
  estimationInit = list(engine = "default_c")
)
print(summary(mod_choice))

# Extract parameters
beta_rate   <- mod_rate$parameters    # intercept, indeg, outdeg
beta_choice <- mod_choice$parameters  # inertia, indeg, outdeg

cat("\nRate parameters:  ", paste(round(beta_rate,   4), collapse = ", "), "\n")
cat("Choice parameters:", paste(round(beta_choice, 4), collapse = ", "), "\n")


# =============================================================================
# 4. EFFECT STATISTICS — standalone functions for simulation
# =============================================================================
# Each function takes:
#   actor_i, actor_j : character labels
#   net_mat          : N×N weight matrix (cumulative, not windowed)
#   recent_events    : data.frame(time, sender, receiver) of events
#                      within the current window
#   current_time     : numeric
# Returns: numeric statistic value

# Helper: windowed event subset
get_window <- function(recent_events, current_time, W) {
  recent_events[recent_events$time > (current_time - W), ]
}

# --- CHOICE effects ---

# inertia(window=W): has i→j occurred in the last W time units?
# s(i,j,t,x) = I(x(t)_ij > 0) in windowed network
stat_inertia_w <- function(i, j, recent, current_time, W) {
  rw <- get_window(recent, current_time, W)
  as.integer(any(rw$sender == i & rw$receiver == j))
}

# indeg(window=W) for choice: how many distinct senders contacted j in window?
# s(i,j,t,x) = sum_k I(x(t)_kj > 0) [windowed, unweighted]
stat_indeg_choice_w <- function(j, recent, current_time, W) {
  rw <- get_window(recent, current_time, W)
  length(unique(rw$sender[rw$receiver == j]))
}

# outdeg(window=W) for choice: how many distinct receivers did j contact in window?
# s(i,j,t,x) = sum_l I(x(t)_jl > 0) [windowed, unweighted, alter type]
stat_outdeg_choice_w <- function(j, recent, current_time, W) {
  rw <- get_window(recent, current_time, W)
  length(unique(rw$receiver[rw$sender == j]))
}

# --- RATE effects ---

# indeg(window=W) for rate: how many senders contacted i in window? [ego type]
stat_indeg_rate_w <- function(i, recent, current_time, W) {
  rw <- get_window(recent, current_time, W)
  length(unique(rw$sender[rw$receiver == i]))
}

# outdeg(window=W) for rate: how many receivers did i contact in window? [ego type]
stat_outdeg_rate_w <- function(i, recent, current_time, W) {
  rw <- get_window(recent, current_time, W)
  length(unique(rw$receiver[rw$sender == i]))
}

# Compute rate for actor i at current time
compute_rate <- function(i, recent, current_time, W, beta_rate) {
  # beta_rate = c(intercept, b_indeg, b_outdeg)
  s_indeg  <- stat_indeg_rate_w(i, recent, current_time, W)
  s_outdeg <- stat_outdeg_rate_w(i, recent, current_time, W)
  linear   <- beta_rate[1] + beta_rate[2] * s_indeg + beta_rate[3] * s_outdeg
  exp(linear)
}

# Compute choice probabilities for sender i over all possible receivers
compute_choice_probs <- function(i, actors, recent, current_time, W, beta_choice) {
  receivers <- actors[actors != i]
  log_utils <- sapply(receivers, function(j) {
    s1 <- stat_inertia_w(i, j, recent, current_time, W)
    s2 <- stat_indeg_choice_w(j, recent, current_time, W)
    s3 <- stat_outdeg_choice_w(j, recent, current_time, W)
    beta_choice[1] * s1 + beta_choice[2] * s2 + beta_choice[3] * s3
  })
  # Numerically stable softmax
  log_utils <- log_utils - max(log_utils)
  probs     <- exp(log_utils)
  probs     / sum(probs)
}


# =============================================================================
# 5. FORWARD SIMULATOR
# =============================================================================
# Simulates one event sequence from the fitted DyNAM.
# Runs until T_max (same time horizon as observed data).
# Returns a data.frame with columns: time, sender, receiver

simulate_dynam <- function(actors, T_max, W,
                           beta_rate, beta_choice,
                           max_events = 10000) {

  sim_events <- data.frame(
    time     = numeric(0),
    sender   = character(0),
    receiver = character(0),
    stringsAsFactors = FALSE
  )

  current_time <- 0

  for (step in seq_len(max_events)) {

    # Draw waiting times for each actor from Exp(rate_i)
    rates <- sapply(actors, function(i) {
      compute_rate(i, sim_events, current_time, W, beta_rate)
    })
    rates <- pmax(rates, 1e-10)  # floor to avoid Inf draws

    waiting_times <- rexp(length(actors), rate = rates)
    names(waiting_times) <- actors

    # Next event: actor with shortest waiting time
    next_actor <- names(which.min(waiting_times))
    next_time  <- current_time + min(waiting_times)

    # Stop if we exceed the observation window
    if (next_time > T_max) break

    # Draw receiver from multinomial choice distribution
    probs    <- compute_choice_probs(next_actor, actors,
                                     sim_events, next_time, W, beta_choice)
    receiver <- sample(actors[actors != next_actor], size = 1, prob = probs)

    # Record event
    sim_events <- rbind(sim_events, data.frame(
      time     = next_time,
      sender   = next_actor,
      receiver = receiver,
      stringsAsFactors = FALSE
    ))

    current_time <- next_time
  }

  sim_events
}


# =============================================================================
# 6. RUN 10 SIMULATIONS
# =============================================================================

n_sim <- 10
cat(sprintf("\nRunning %d simulations...\n", n_sim))

sim_list <- vector("list", n_sim)
pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

for (s in seq_len(n_sim)) {
  sim_list[[s]] <- simulate_dynam(
    actors      = all_actors,
    T_max       = T_max,
    W           = W,
    beta_rate   = beta_rate,
    beta_choice = beta_choice
  )
  setTxtProgressBar(pb, s)
}
close(pb)

n_events_sim <- sapply(sim_list, nrow)
cat(sprintf("\nSimulated events: mean=%.1f, sd=%.1f, obs=%d\n",
            mean(n_events_sim), sd(n_events_sim), nrow(events_df)))


# =============================================================================
# 7. NETWORK SUMMARY STATISTICS — for GoF
# =============================================================================

# From an event list, compute summary statistics for GoF
network_stats <- function(ev, actors) {

  if (nrow(ev) == 0) {
    return(list(
      indeg       = setNames(rep(0, length(actors)), actors),
      outdeg      = setNames(rep(0, length(actors)), actors),
      tie_weights = numeric(0),
      recip_ratio = 0,
      w_density   = 0
    ))
  }

  # Indegree and outdegree (number of distinct partners)
  indeg <- table(factor(ev$receiver, levels = actors))
  outdeg <- table(factor(ev$sender,  levels = actors))

  # Tie weight per dyad (total events per directed dyad)
  tie_w <- ev %>%
    count(sender, receiver, name = "weight")

  # Reciprocity: fraction of directed dyads with a return tie
  # For each (i,j) with weight > 0, check if (j,i) also exists
  dyads <- paste(tie_w$sender, tie_w$receiver, sep = "->")
  rev_dyads <- paste(tie_w$receiver, tie_w$sender, sep = "->")
  recip_flag <- rev_dyads %in% dyads
  recip_ratio <- mean(recip_flag)

  # Weighted density: total events / (N*(N-1)) possible directed pairs
  N  <- length(actors)
  w_density <- nrow(ev) / (N * (N - 1))

  list(
    indeg       = as.numeric(indeg),
    outdeg      = as.numeric(outdeg),
    tie_weights = tie_w$weight,
    recip_ratio = recip_ratio,
    w_density   = w_density
  )
}

# Observed statistics
obs_stats <- network_stats(events_df, all_actors)

# Simulated statistics
sim_stats <- lapply(sim_list, network_stats, actors = all_actors)


# =============================================================================
# 8. GOF PLOTS
# =============================================================================

# Helper: extract one statistic across all simulations
extract_sim_stat <- function(sim_stats, stat_name) {
  lapply(sim_stats, `[[`, stat_name)
}

# --- 8A. Indegree distribution ---
sim_indeg <- extract_sim_stat(sim_stats, "indeg")
obs_indeg <- obs_stats$indeg

# Build distribution comparison
max_deg <- max(c(unlist(sim_indeg), obs_indeg))
deg_vals <- 0:max_deg

sim_indeg_dist <- do.call(rbind, lapply(seq_along(sim_indeg), function(s) {
  counts <- table(factor(sim_indeg[[s]], levels = deg_vals))
  data.frame(sim = s, degree = deg_vals, count = as.numeric(counts))
}))

obs_indeg_dist <- data.frame(
  degree = deg_vals,
  count  = as.numeric(table(factor(obs_indeg, levels = deg_vals)))
)

# Summarize simulated distribution
sim_indeg_summary <- sim_indeg_dist %>%
  group_by(degree) %>%
  summarise(
    median = median(count),
    q05    = quantile(count, 0.05),
    q95    = quantile(count, 0.95),
    .groups = "drop"
  )

p_indeg <- ggplot() +
  geom_ribbon(data = sim_indeg_summary,
              aes(x = degree, ymin = q05, ymax = q95),
              fill = "#3A7DC9", alpha = 0.25) +
  geom_line(data = sim_indeg_summary,
            aes(x = degree, y = median),
            colour = "#3A7DC9", linewidth = 0.8) +
  geom_point(data = obs_indeg_dist,
             aes(x = degree, y = count),
             colour = "#E8735A", size = 2.5) +
  geom_line(data = obs_indeg_dist,
            aes(x = degree, y = count),
            colour = "#E8735A", linewidth = 0.8) +
  labs(
    title    = "GoF: Indegree distribution",
    subtitle = "Blue = simulated (median + 90% CI) | Red = observed",
    x        = "Indegree", y = "Count"
  ) +
  theme_minimal(base_size = 12)
print(p_indeg)

# --- 8B. Outdegree distribution ---
sim_outdeg <- extract_sim_stat(sim_stats, "outdeg")
obs_outdeg <- obs_stats$outdeg

max_deg2 <- max(c(unlist(sim_outdeg), obs_outdeg))
deg_vals2 <- 0:max_deg2

sim_outdeg_dist <- do.call(rbind, lapply(seq_along(sim_outdeg), function(s) {
  counts <- table(factor(sim_outdeg[[s]], levels = deg_vals2))
  data.frame(sim = s, degree = deg_vals2, count = as.numeric(counts))
}))

obs_outdeg_dist <- data.frame(
  degree = deg_vals2,
  count  = as.numeric(table(factor(obs_outdeg, levels = deg_vals2)))
)

sim_outdeg_summary <- sim_outdeg_dist %>%
  group_by(degree) %>%
  summarise(
    median = median(count),
    q05    = quantile(count, 0.05),
    q95    = quantile(count, 0.95),
    .groups = "drop"
  )

p_outdeg <- ggplot() +
  geom_ribbon(data = sim_outdeg_summary,
              aes(x = degree, ymin = q05, ymax = q95),
              fill = "#3A7DC9", alpha = 0.25) +
  geom_line(data = sim_outdeg_summary,
            aes(x = degree, y = median),
            colour = "#3A7DC9", linewidth = 0.8) +
  geom_point(data = obs_outdeg_dist,
             aes(x = degree, y = count),
             colour = "#E8735A", size = 2.5) +
  geom_line(data = obs_outdeg_dist,
            aes(x = degree, y = count),
            colour = "#E8735A", linewidth = 0.8) +
  labs(
    title    = "GoF: Outdegree distribution",
    subtitle = "Blue = simulated (median + 90% CI) | Red = observed",
    x        = "Outdegree", y = "Count"
  ) +
  theme_minimal(base_size = 12)
print(p_outdeg)

# --- 8C. Tie weight distribution ---
sim_tw <- extract_sim_stat(sim_stats, "tie_weights")
obs_tw <- obs_stats$tie_weights

max_tw <- max(c(unlist(sim_tw), obs_tw), na.rm = TRUE)
tw_vals <- 1:max_tw

sim_tw_dist <- do.call(rbind, lapply(seq_along(sim_tw), function(s) {
  tw <- sim_tw[[s]]
  if (length(tw) == 0) tw <- integer(0)
  counts <- table(factor(tw, levels = tw_vals))
  data.frame(sim = s, weight = tw_vals, count = as.numeric(counts))
}))

obs_tw_dist <- data.frame(
  weight = tw_vals,
  count  = as.numeric(table(factor(obs_tw, levels = tw_vals)))
)

sim_tw_summary <- sim_tw_dist %>%
  group_by(weight) %>%
  summarise(
    median = median(count),
    q05    = quantile(count, 0.05),
    q95    = quantile(count, 0.95),
    .groups = "drop"
  )

p_tw <- ggplot() +
  geom_ribbon(data = sim_tw_summary,
              aes(x = weight, ymin = q05, ymax = q95),
              fill = "#3A7DC9", alpha = 0.25) +
  geom_line(data = sim_tw_summary,
            aes(x = weight, y = median),
            colour = "#3A7DC9", linewidth = 0.8) +
  geom_point(data = obs_tw_dist,
             aes(x = weight, y = count),
             colour = "#E8735A", size = 2.5) +
  geom_line(data = obs_tw_dist,
            aes(x = weight, y = count),
            colour = "#E8735A", linewidth = 0.8) +
  labs(
    title    = "GoF: Tie weight distribution (per dyad)",
    subtitle = "Blue = simulated (median + 90% CI) | Red = observed",
    x        = "Tie weight (total events per dyad)", y = "Count"
  ) +
  theme_minimal(base_size = 12)
print(p_tw)

# --- 8D. Reciprocity distribution ---
sim_recip <- sapply(sim_stats, `[[`, "recip_ratio")
obs_recip <- obs_stats$recip_ratio

p_recip <- ggplot(data.frame(recip = sim_recip), aes(x = recip)) +
  geom_histogram(bins = 30, fill = "#3A7DC9", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = obs_recip, colour = "#E8735A",
             linewidth = 1.2, linetype = "solid") +
  annotate("text", x = obs_recip, y = Inf, vjust = 2, hjust = -0.1,
           label = sprintf("Observed\n%.3f", obs_recip),
           colour = "#E8735A", size = 3.5) +
  labs(
    title    = "GoF: Reciprocity (fraction of reciprocated dyads)",
    subtitle = "Blue histogram = simulated | Red line = observed",
    x        = "Reciprocity ratio", y = "Count"
  ) +
  theme_minimal(base_size = 12)
print(p_recip)

# --- 8E. Weighted density ---
sim_wdens <- sapply(sim_stats, `[[`, "w_density")
obs_wdens <- obs_stats$w_density

p_wdens <- ggplot(data.frame(wd = sim_wdens), aes(x = wd)) +
  geom_histogram(bins = 30, fill = "#3A7DC9", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = obs_wdens, colour = "#E8735A",
             linewidth = 1.2) +
  annotate("text", x = obs_wdens, y = Inf, vjust = 2, hjust = -0.1,
           label = sprintf("Observed\n%.3f", obs_wdens),
           colour = "#E8735A", size = 3.5) +
  labs(
    title    = "GoF: Weighted density (total events / N*(N-1))",
    subtitle = "Blue histogram = simulated | Red line = observed",
    x        = "Weighted density", y = "Count"
  ) +
  theme_minimal(base_size = 12)
print(p_wdens)


# =============================================================================
# 9. SUMMARY TABLE
# =============================================================================

gof_summary <- data.frame(
  Statistic   = c("Reciprocity", "Weighted density",
                  "Mean indegree", "SD indegree",
                  "Mean outdegree", "SD outdegree",
                  "Mean tie weight", "SD tie weight",
                  "N events"),
  Observed    = c(
    obs_stats$recip_ratio,
    obs_stats$w_density,
    mean(obs_stats$indeg),
    sd(obs_stats$indeg),
    mean(obs_stats$outdeg),
    sd(obs_stats$outdeg),
    ifelse(length(obs_stats$tie_weights) > 0,
           mean(obs_stats$tie_weights), NA),
    ifelse(length(obs_stats$tie_weights) > 0,
           sd(obs_stats$tie_weights), NA),
    nrow(events_df)
  ),
  Sim_mean    = c(
    mean(sapply(sim_stats, `[[`, "recip_ratio")),
    mean(sapply(sim_stats, `[[`, "w_density")),
    mean(sapply(sim_stats, function(s) mean(s$indeg))),
    mean(sapply(sim_stats, function(s) sd(s$indeg))),
    mean(sapply(sim_stats, function(s) mean(s$outdeg))),
    mean(sapply(sim_stats, function(s) sd(s$outdeg))),
    mean(sapply(sim_stats, function(s)
      ifelse(length(s$tie_weights)>0, mean(s$tie_weights), NA)), na.rm=TRUE),
    mean(sapply(sim_stats, function(s)
      ifelse(length(s$tie_weights)>0, sd(s$tie_weights), NA)), na.rm=TRUE),
    mean(n_events_sim)
  ),
  Sim_sd      = c(
    sd(sapply(sim_stats, `[[`, "recip_ratio")),
    sd(sapply(sim_stats, `[[`, "w_density")),
    sd(sapply(sim_stats, function(s) mean(s$indeg))),
    sd(sapply(sim_stats, function(s) sd(s$indeg))),
    sd(sapply(sim_stats, function(s) mean(s$outdeg))),
    sd(sapply(sim_stats, function(s) sd(s$outdeg))),
    sd(sapply(sim_stats, function(s)
      ifelse(length(s$tie_weights)>0, mean(s$tie_weights), NA)), na.rm=TRUE),
    sd(sapply(sim_stats, function(s)
      ifelse(length(s$tie_weights)>0, sd(s$tie_weights), NA)), na.rm=TRUE),
    sd(n_events_sim)
  )
)

cat("\n=== GoF Summary ===\n")
print(gof_summary, digits = 3, row.names = FALSE)


# =============================================================================
# 10. SAVE
# =============================================================================

ggsave(file.path(out_dir, "gof_indegree.png"),   p_indeg,  width=7, height=5, dpi=300)
ggsave(file.path(out_dir, "gof_outdegree.png"),  p_outdeg, width=7, height=5, dpi=300)
ggsave(file.path(out_dir, "gof_tieweights.png"), p_tw,     width=7, height=5, dpi=300)
ggsave(file.path(out_dir, "gof_reciprocity.png"),p_recip,  width=7, height=5, dpi=300)
ggsave(file.path(out_dir, "gof_wdensity.png"),   p_wdens,  width=7, height=5, dpi=300)

save(mod_rate, mod_choice, sim_list, sim_stats, obs_stats, gof_summary,
     file = file.path(out_dir, "dynam_gof_flies.RData"))

cat("\nAll outputs saved to:", out_dir, "\n")
