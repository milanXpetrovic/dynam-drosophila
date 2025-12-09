# Package Installation
# install.packages(c('here', 'remotes', 'glue', 'cmdstanr', 'foreach', 'doParallel'))
# remotes::install_github("stocnet/goldfish@develop")
# remotes::install_github('snlab-ch/goldfish.latent@develop')

# Package Loading
library(goldfish.latent)
library(goldfish)
library(here)
library(glue)
library(cmdstanr)
library(foreach)
library(doParallel)
library(parallel)
library(iterators)

# Utility Functions
source(here::here("src/utils.R"))

# Global Parameters
treatments_list <- c('CS_10D', 'Cs_5DIZ', 'CsCh')

base_path <- "/home/milky/drosophila-pipelnie/6_dynam_models"
edgelist_path <- file.path(base_path, "data", "edgelists")
covariances_path <- file.path(base_path, "data", "covariances")
soc_space_matrix_path <- file.path(base_path, "data", "soc_space_matrix")
results_path <- file.path(base_path, "res")
results_rate_model_path <- file.path(results_path, "rate_model_individual")
results_choice_model_path <- file.path(results_path, "choice_model_individual")

dir.create(results_rate_model_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_choice_model_path, showWarnings = FALSE, recursive = TRUE)

choice_model_formula <- dependent ~
    # indeg(interaction_network, weighted = TRUE, window = 96) +
    # indeg(interaction_network, weighted = TRUE, window = 288) +
    # indeg(interaction_network, weighted = TRUE, window = 864) +
    # indeg(interaction_network, weighted = TRUE, window = 2592) +

    # outdeg(interaction_network, weighted = TRUE, window = 96) +
    # outdeg(interaction_network, weighted = TRUE, window = 288) +
    # outdeg(interaction_network, weighted = TRUE, window = 864) +
    # outdeg(interaction_network, weighted = TRUE, window = 2592) +

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

    cycle(interaction_network, window = 96) +
    cycle(interaction_network, window = 288) +
    cycle(interaction_network, window = 864) +
    cycle(interaction_network, window = 2592) +
    
    inertia(interaction_network, weighted = TRUE, window = 96) +
    inertia(interaction_network, weighted = TRUE, window = 288) +
    inertia(interaction_network, weighted = TRUE, window = 864) +
    inertia(interaction_network, weighted = TRUE, window = 2592) +

    recip(interaction_network, weighted = TRUE, window = 96) +
    recip(interaction_network, weighted = TRUE, window = 288) +
    recip(interaction_network, weighted = TRUE, window = 864) +
    recip(interaction_network, weighted = TRUE, window = 2592) +

    inertia(interaction_network, weighted = FALSE, window = 96) +
    # inertia(interaction_network, weighted = FALSE, window = 288) +
    inertia(interaction_network, weighted = FALSE, window = 864) +
    inertia(interaction_network, weighted = FALSE, window = 2592) +

    recip(interaction_network, weighted = FALSE, window = 96) +
    recip(interaction_network, weighted = FALSE, window = 288) +
    recip(interaction_network, weighted = FALSE, window = 864) +
    recip(interaction_network, weighted = FALSE, window = 2592) +

    common_sender(interaction_network, window = 96) +
    common_sender(interaction_network, window = 288) +
    common_sender(interaction_network, window = 864) +
    common_sender(interaction_network, window = 2592) +

    common_receiver(interaction_network, window = 96) +
    common_receiver(interaction_network, window = 288) +
    common_receiver(interaction_network, window = 864) +
    common_receiver(interaction_network, window = 2592) +

    alter(nodesAttr$positive_influence) +
    alter(nodesAttr$negative_influence) +
    alter(nodesAttr$positive_inf_weighted) +
    alter(nodesAttr$negative_inf_weighted) +

    alter(nodesAttr$distance_traveled_between_interactions) + 
    alter(nodesAttr$number_of_flies_in_soc_space) + 
    alter(nodesAttr$unique_partners_met_interaction_space) + 
    alter(nodesAttr$unique_partners_met_social_space) + 

    sim(nodesAttr$activity) +
    sim(nodesAttr$popularity) + 

    alter(flies$young) +
    alter(flies$old) +
    alter(flies$isolated)

log_file_path <- file.path(results_path, "parallel_execution_log.txt")
message(paste("Parallel execution log will be written to:", log_file_path))

all_files_list <- lapply(treatments_list, function(treatment) {
    folder_path <- file.path(edgelist_path, treatment)
    csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
    if (length(csv_files) > 0) {
        data.frame(
            treatment_name = treatment,
            file_path = csv_files,
            group_name = tools::file_path_sans_ext(basename(csv_files))
        )
    } else {
        NULL
    }
})

files_to_process <- do.call(rbind, all_files_list)
# files_to_process <- head(files_to_process, 4) ##! TAKE CARE 

n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

parallel_results <- foreach(
    i = 1:nrow(files_to_process),
    .packages = c("goldfish", "goldfish.latent", "glue", "tools"),
    # .export = c("read_and_clean", "choice_model_formula", "covariances_path", "soc_space_matrix_path" ),
    .errorhandling = "pass"
) %dopar% {
    worker_pid <- Sys.getpid()
    start_time <- Sys.time()

    start_message <- glue::glue(
        "{format(start_time, '%Y-%m-%d %H:%M:%S')} [PID: {worker_pid}] START - Processing group: {files_to_process$group_name[i]}"
    )
    
    cat(paste0(start_message, "\n"), file = log_file_path, append = TRUE)

    # Per-group variables
    # i <- 1
    treatment_name <- files_to_process$treatment_name[i]
    group_name <- files_to_process$group_name[i]
    file_path <- files_to_process$file_path[i]
    # Load interaction data
    interaction_data <- read.csv(file_path)
    interaction_data <- interaction_data[, c('time', 'sender', 'receiver', 'increment')] 
    # interaction_data <- head(interaction_data)
    # Initialize flies data frame with attributes
    flies <- getActors(interaction_data)
    flies$popularity <- 0
    flies$activity <- 0
    flies$popularity_weighted <- 0
    flies$activity_weighted <- 0
    flies$positive_influence <- 0
    flies$negative_influence <- 0
    flies$positive_inf_weighted <- 0
    flies$negative_inf_weighted <- 0
    flies$distance_traveled_between_interactions <- 0
    flies$number_of_flies_in_soc_space <- 0
    flies$unique_partners_met_interaction_space <- 0
    flies$unique_partners_met_social_space <- 0
    flies$young <- 0
    flies$old <- 0
    flies$isolated <- 0

    if (treatment_name == "CsCh") flies$young <- 1
    if (treatment_name == "CS_10D") flies$old <- 1
    if (treatment_name == "Cs_5DIZ") flies$isolated <- 1

    # Load covariance data
    group_cov_path <- file.path(covariances_path, treatment_name, group_name)
    cov_in_degree <- read_and_clean('in_degree.csv', group_cov_path)
    cov_out_degree <- read_and_clean('out_degree.csv', group_cov_path)
    cov_in_weighted <- read_and_clean('in_weighted.csv', group_cov_path)
    cov_out_weighted <- read_and_clean('out_weighted.csv', group_cov_path)
    cov_positiveinfluence <- read_and_clean('positiveinfluence.csv', group_cov_path)
    cov_negativeinfluence <- read_and_clean('negativeinfluence.csv', group_cov_path)
    cov_positiveinf_weighted <- read_and_clean('positiveinf_weighted.csv', group_cov_path)
    cov_negativeinf_weighted <- read_and_clean('negativeinf_weighted.csv', group_cov_path)
    cov_distance_traveled_between_interactions <- read_and_clean('distance_traveled_between_interactions.csv', group_cov_path)
    cov_number_of_flies_in_soc_space <- read_and_clean('number_of_flies_in_soc_space.csv', group_cov_path)
    cov_unique_partners_met_interaction_space <- read_and_clean('unique_partners_met_interaction_space.csv', group_cov_path)
    cov_unique_partners_met_social_space <- read_and_clean('unique_partners_met_social_space.csv', group_cov_path)
    
    # Create and link goldfish objects
    nodesAttr <- make_nodes(flies)
    interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
    interaction_network <- link_events(x = interaction_network, change_event = interaction_data, nodes = nodesAttr)
    dependent <- make_dependent_events(events = interaction_data, nodes = nodesAttr, default_network = interaction_network)

    # Link covariance attributes
    nodesAttr <- link_events(nodesAttr, cov_in_degree, 'popularity')
    nodesAttr <- link_events(nodesAttr, cov_out_degree, 'activity')
    nodesAttr <- link_events(nodesAttr, cov_in_weighted, 'popularity_weighted')
    nodesAttr <- link_events(nodesAttr, cov_out_weighted, 'activity_weighted')
    nodesAttr <- link_events(nodesAttr, cov_positiveinfluence, 'positive_influence')
    nodesAttr <- link_events(nodesAttr, cov_negativeinfluence, 'negative_influence')
    nodesAttr <- link_events(nodesAttr, cov_positiveinf_weighted, 'positive_inf_weighted')
    nodesAttr <- link_events(nodesAttr, cov_negativeinf_weighted, 'negative_inf_weighted')
    nodesAttr <- link_events(nodesAttr, cov_distance_traveled_between_interactions, 'distance_traveled_between_interactions')
    nodesAttr <- link_events(nodesAttr, cov_number_of_flies_in_soc_space, 'number_of_flies_in_soc_space')
    nodesAttr <- link_events(nodesAttr, cov_unique_partners_met_interaction_space, 'unique_partners_met_interaction_space')
    nodesAttr <- link_events(nodesAttr, cov_unique_partners_met_social_space, 'unique_partners_met_social_space')

    # # Load and process social space matrix
    # group_ssm_path <- file.path(soc_space_matrix_path, treatment_name, group_name)
    # initial_social_space_matrix <- as.matrix(read.csv(file.path(group_ssm_path, "m_at_t0.csv"), row.names = 1))
    # social_space_updates <- read.csv(file.path(group_ssm_path, "updates.csv"))
    # if (colnames(social_space_updates)[1] == "X") social_space_updates <- social_space_updates[, -1]
    # names(social_space_updates) <- c("time", "sender", "receiver", "replace")
    # social_space_updates$sender <- glue("fly{sender+1}", sender = social_space_updates$sender)
    # social_space_updates$receiver <- glue("fly{receiver+1}", receiver = social_space_updates$receiver)
    
    # # Create and link social space network
    # initial_social_space_matrix <- initial_social_space_matrix[flies$label, flies$label]
    # socialSpace_network <- make_network(initial_social_space_matrix, nodes = nodesAttr, directed = FALSE)
    # socialSpace_network <- link_events(socialSpace_network, social_space_updates, nodes = nodesAttr)

    # Prepare data for modeling
    dataDynam <- make_data(nodesAttr, interaction_network, dependent, flies)
    
    print(glue("Processing group: {group_name}"))
    data_for_this_group_choice <- make_data_re(
        random_effects = list(inertia(interaction_network, weighted = FALSE, window = 288) ~ 1),
        fixed_effects = choice_model_formula,
        # support_constraint = ~ tie(socialSpace_network),
        model = "DyNAM",
        sub_model = "choice",
        data = dataDynam,
        # control_preprocessing = set_preprocessing_opt(start_time = 0, end_time = 2020)
    )

    # Sample data
    samplesGroups <- sample_data(data_for_this_group_choice, sample_size_set = 1)
    end_time <- Sys.time()
    duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    end_message <- glue::glue(
        "{format(end_time, '%Y-%m-%d %H:%M:%S')} [PID: {worker_pid}] END   - Finished group: {group_name}. Duration: {duration}s"
    )
    cat(paste0(end_message, "\n"), file = log_file_path, append = TRUE)

    # Return results for this group
    list(
        group_name = group_name,
        treatment_name = treatment_name,
        group_index = i,
        data_object = data_for_this_group_choice,
        samples_object = samplesGroups
    )
}

# Stop the parallel cluster
stopCluster(cl)
registerDoSEQ()

all_groups_data_list <- list()
all_groups_samples <- list()
all_groups_info_list <- list()

for (result in parallel_results) {
    if (inherits(result, "list")) { # Check if the result is a list (not an error)
        group_name <- result$group_name
        all_groups_data_list[[group_name]] <- result$data_object
        all_groups_samples[[group_name]] <- result$samples_object
        all_groups_info_list[[group_name]] <- data.frame(
            group = group_name,
            treatment = result$treatment_name,
            ixGroup = result$group_index
        )
    } else {
        # Log or print the error from the parallel task
        print(paste("An error occurred for a group:", result))
    }
}

all_groups_info <- do.call(rbind, all_groups_info_list)

saveRDS(all_groups_data_list, file = "all_groups_data_list.rds")
saveRDS(all_groups_samples, file = "all_groups_samples.rds")
saveRDS(all_groups_info, file = "all_groups_info.rds")
print("Processing complete. Results saved.")

########################################################################################################################

all_groups_data_list <- readRDS("tmp/all_groups_data_list.rds")
all_groups_info <- readRDS("tmp/all_groups_info.rds")

CORRELATION_THRESHOLD <- 0.75

all_colnames <- unique(unlist(lapply(all_groups_data_list, function(group) {
    colnames(group$data_stan$X_choice[, -1, drop = FALSE])})))
all_colnames <- sort(all_colnames)

final_count_matrix <- matrix(0,
    nrow = length(all_colnames),
    ncol = length(all_colnames),
    dimnames = list(all_colnames, all_colnames)
)

for (group_name in names(all_groups_data_list)) {
    fixed_effects_matrix <- all_groups_data_list[[group_name]]$data_stan$X_choice[, -1, drop = FALSE]
    constant_cols <- which(apply(fixed_effects_matrix, 2, var) == 0)
    if (length(constant_cols) > 0) {
        fixed_effects_matrix <- fixed_effects_matrix[, -constant_cols, drop = FALSE]
    }
    if (ncol(fixed_effects_matrix) < 2) {
        next
    }
    cor_matrix_group <- cor(fixed_effects_matrix)
    high_corr_adj_matrix <- (abs(cor_matrix_group) > CORRELATION_THRESHOLD) * 1
    diag(high_corr_adj_matrix) <- 0
    current_predictors <- colnames(high_corr_adj_matrix)
    final_count_matrix[current_predictors, current_predictors] <-
        final_count_matrix[current_predictors, current_predictors] + high_corr_adj_matrix
}

print(final_count_matrix)

# install.packages("corrplot")
library(corrplot)
library(ggplot2)
library(reshape2)

color_palette <- colorRampPalette(c("white", "red"))(100)

relevant_vars_idx <- (rowSums(final_count_matrix) > 0) | (colSums(final_count_matrix) > 0)
matrix_to_plot <- final_count_matrix[relevant_vars_idx, relevant_vars_idx]
png("high_correlation_plot.png", width = 14, height = 14, units = "in", res = 300)
corrplot(
    matrix_to_plot,
    method = "color",
    col = color_palette,
    type = "upper",
    addgrid.col = "grey",
    is.corr = FALSE,
    cl.lim = c(0, 60),
    tl.col = "black",
    tl.srt = 90,        # Rotate labels
    tl.cex = 0.9,       # Adjust label font size
    addCoef.col = "black",
    number.cex = 0.01,
    title = CORRELATION_THRESHOLD,
    mar = c(0, 0, 1, 0)
)
dev.off()

############################################################################################################

source(here::here("src/utils.R"))
data_stan_re <- GatherGroups(
    groupsInfo = all_groups_info,
    data = all_groups_data_list
)

attr(data_stan_re, "sample") <- FALSE
dataCox <- goldfish.latent:::make_cox_model_data(data_stan_re, "choice")
names(dataCox) <- stringr::str_replace_all(names(dataCox), ":", "__")
termsAdd <- names(dataCox)[seq.int(which(names(dataCox) == "event") - 1)]

# cox model
##! full dataset for the group
mod00 <- coxph(
    as.formula(paste0("Surv(uno, selected) ~", paste(termsAdd, collapse = " + "),
        "+ strata(event)")),
    data = dataCox,
    control = coxph.control(iter.max = 1e3)
    )

library(broom)
coefTable <- tidy(mod00, conf.int = TRUE) 
write.csv(coefTable, file = "res.csv")

termsInteract <- names(dataCox)[seq.int(which(names(dataCox) == "sim_nodesAttrOfpopularity") - 1)]
termsExpCond <- grep("alter_fliesOf", names(dataCox), value = TRUE)
interaction_formula <- as.formula(paste0(
    "Surv(uno, selected) ~ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ") + strata(event)"
))
mod01 <- coxph(interaction_formula,
    data = dataCox,
    control = coxph.control(iter.max = 1e3)
    )
coefTable01 <- tidy(mod01, conf.int = TRUE)
write.csv(coefTable01, file = "res01.csv")


dataCox[, grep("distance", names(dataCox), value = T)] <- log1p(dataCox[, grep("distance", names(dataCox), value = T)])
termsInteract <- names(dataCox)[seq.int(which(names(dataCox) == "sim_nodesAttrOfpopularity") - 1)]
termsInteract <- termsInteract[!grepl("W_96$", termsInteract)]
termsExpCond <- grep("alter_fliesOf", names(dataCox), value = TRUE)
interaction_formula <- as.formula(paste0(
    "Surv(uno, selected) ~ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ") + strata(event)"
))
mod02 <- coxph(interaction_formula,
    data = dataCox,
    control = coxph.control(iter.max = 1e3)
    )
coefTable02 <- tidy(mod02, conf.int = TRUE)
write.csv(coefTable02, file = "tmp/res02.csv")


all_groups_cox_coeff_list[[group_name]] <- tidy(mod00)

##! sample
dataCox <- make_cox_model_data(samplesGroups, "choice")
names(dataCox) <- stringr::str_replace_all(names(dataCox), ":", "__")
termsAdd <- names(dataCox)[seq.int(which(names(dataCox) == "event") - 1)]

mod00 <- coxph(
    as.formula(paste0("Surv(uno, selected) ~", paste(termsAdd, collapse = " + "),
        "+ strata(event) + offset(offset_sample)")),
    data = dataCox,
    control = coxph.control(iter.max = 1e3)
    )



all_groups_sample_cox_coeff_list[[group_name]] <- tidy(mod00)

### stan bro
# install.packages("tidy")
library(survival)
# library(tidy)

all_groups_data_list_fixed <- lapply(all_groups_data_list, function(group_data) {
  names(group_data)[names(group_data) == "data_stan"] <- "dataStan"
  return(group_data)
})
##! Problems with naming data_stan and dataStan

source(here::here("src/utils.R"))
data_stan_re <- GatherGroups(
    groupsInfo = all_groups_info,
    data = all_groups_data_list
)

attr(data_stan_re, "subModel") <- "choice"
stanCode <- make_model_code(data_stan_re)
modRandomEffects <- cmdstan_model(stanCode)

random_effects_model <- modRandomEffects$sample(
    data = data_stan_re$data_stan,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500
)

print(summary(random_effects_model))
model_save_path <- file.path(results_path, "full_random_effects_model.rds")
saveRDS(random_effects_model, file = model_save_path)

modRandomEffects <- cmdstan_model('./src/dynam_map_reduce_event.stan', cpp_options = list(stan_threads = TRUE), force_recompile=TRUE)

data_stan_re$data_stan$grain_size <- 4

old_X_choice <- data_stan_re$data_stan$X_choice
namesVars <- colnames(data_stan_re$data_stan$X_choice)
termsInteract <- namesVars[seq.int(which(namesVars == "sim_nodesAttrOfpopularity") - 1)]
termsExpCond <- grep("alter_fliesOf", namesVars, value = TRUE)
int_formula <- as.formula(paste0(
    " ~ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ") + 0"
))

var_distance <- grep("distance", colnames(old_X_choice), value = T)
old_X_choice[, var_distance] <- log1p(old_X_choice[, var_distance])
new_X_choice <- model.matrix(int_formula, data = as.data.frame(old_X_choice))
new_X_scale <- scale(new_X_choice)
data_stan_re$scale <- list( 
    center = attr(new_X_scale, "scaled:center"),
    scale = attr(new_X_scale, "scaled:scale")
)

Z_scale <- scale(data_stan_re$data_stan$Z_choice)
data_stan_re$scale_Z <- list( 
    center = attr(Z_scale, "scaled:center"),
    scale = attr(Z_scale, "scaled:scale")
)

# data_stan_re$data_stan$Z_choice <- data_stan_re$data_stan$X_choice[, "inertia_interaction_network_2592", drop = FALSE]
data_stan_re$data_stan$Z_choice <- Z_scale[, , drop = FALSE]
data_stan_re$data_stan$X_choice <- new_X_choice[, ]
data_stan_re$data_stan$P_choice <- ncol(new_X_choice)

#!  73807.6 seconds. / 1230 minutes / cca 21h
random_effects_model <- modRandomEffects$sample(
    data = data_stan_re$data_stan,
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 11,
    iter_warmup = 500,
    iter_sampling = 500
)

random_effects_model$save_object(file = "tmp/full_random_effects_model.rds")
arrayDraws <- random_effects_model$draws()
saveRDS(arrayDraws, file = "tmp/arrayDraws.rds")

#### without weighted short window
data_stan_wsw <- data_stan_re$data_stan
data_stan_wsw$X_choice <- new_X_choice[, !grepl("W_96$", colnames(new_X_choice))]
data_stan_wsw$P_choice <- ncol(data_stan_wsw$X_choice)
summary(new_X_choice[, grep("distance", colnames(new_X_choice), value = T)])

random_effects_model_wsw <- modRandomEffects$sample(
    data = data_stan_wsw,
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 11,
    iter_warmup = 500,
    iter_sampling = 500
)

print(summary(random_effects_model_wsw))
model_save_path <- file.path(results_path, "tmp/full_random_effects_model_wsw.rds")
saveRDS(random_effects_model_wsw, file = model_save_path)

draws_beta <- random_effects_model_wsw$draws("beta_choice", format = "draws_matrix")

colnames(draws_beta) <- colnames(data_stan_wsw$X_choice)

namesCols = cbind(beta = colnames(draws_beta), origin = colnames(data_stan_wsw$data_stan$X_choice))
namesPlot = namesCols[grepl("(cycle|trans).+96", namesCols[, 2]), ]

beta_intervals_original <- mcmc_intervals_data(draws_beta)
write.csv(beta_intervals_original, file = "tmp/beta_intervals_original_wsw2.csv")

all_groups_info <- readRDS("tmp/all_groups_info.rds")
groups <- all_groups_info$group
write.csv(groups, file = "tmp/groups.csv")

gamma_intervals <- readRDS("tmp/gamma_intervals_wsw.rds")

draws_gamma <- random_effects_model_wsw$draws(c("gamma", "sigma"), format = "draws_matrix")
# draws_gamma <- draws_gamma / data_stan_re$scale_Z$scale

colnames(draws_gamma) <- colnames(data_stan_wsw$X_choice)
all_groups_info

gamma_intervals <- mcmc_intervals_data(draws_gamma)
write.csv(gamma_intervals, file = "tmp/gamma_intervals_wsw.csv")

#!   chain_id  warmup sampling   total
#! 1        1 12976.2  4632.93 17609.1
random_effects_model_wsw$time()
###########################################################################################

random_effects_model$time()
library(posterior)
library(bayesplot)
library(ggplot2)

sumDraws <- random_effects_model$summary()
random_effects_model$diagnostic_summary()

color_scheme_set("brightblue")
rhats <- rhat(random_effects_model)
rhats=mcmc_rhat(rhats)
ggsave("tmp/rhats.png", rhats, width = 7, height = 10)

ratios_ess <- neff_ratio(random_effects_model)
ratios_ess_plot <- mcmc_neff(ratios_ess, size = 2)
ggsave("tmp/ratios_ess.png", ratios_ess_plot, width = 7, height = 10)

random_effects_model$summary("sigma")

trace_some <- mcmc_trace(random_effects_model$draws(), regex_pars = "beta_choice\\W\\d\\W$")
ggsave("tmp/trace_1.png", trace_some, width = 7, height = 10)

trace_random <- mcmc_trace(random_effects_model$draws(), regex_pars = "sigma|gamma\\W\\d\\W$")
ggsave("tmp/trace_2.png", trace_random, width = 7, height = 10)

caterpillarplot <- mcmc_intervals(random_effects_model$draws("gamma"))
ggsave("tmp/catterpillarplot.png", caterpillarplot, width = 7, height = 10)

# rescale
draws_beta <- random_effects_model$draws("beta_choice", format = "draws_matrix")
draws_beta_rescale <- sweep(draws_beta, 2, data_stan_re$scale$scale, FUN = '/' )
colnames(draws_beta) <- colnames(data_stan_re$data_stan$X_choice)

namesCols = cbind(beta = colnames(draws_beta), origin = colnames(data_stan_re$data_stan$X_choice))
namesPlot = namesCols[grepl("(cycle|trans).+96", namesCols[, 2]), ]

plotTrans96 = mcmc_intervals(random_effects_model$draws(namesPlot[, 1]))
ggsave("tmp/intervals_trans_96.png", plotTrans96, width = 7, height = 10)


draws_gamma <- random_effects_model$draws(c("gamma", "sigma"), format = "draws_matrix")
draws_gamma <- draws_gamma / data_stan_re$scale_Z$scale

beta_intervals <- mcmc_intervals_data(draws_beta_rescale)
write.csv(beta_intervals, file = "tmp/beta_intervals.csv")

beta_intervals_original <- mcmc_intervals_data(draws_beta)
write.csv(beta_intervals_original, file = "tmp/beta_intervals_original.csv")

gamma_intervals <- mcmc_intervals_data(draws_gamma)
write.csv(gamma_intervals, file = "tmp/gamma_intervals.csv")


########################################################################################################################
corr_matrix <- cor(data_stan_re$data_stan$X_choice[, grepl("inertia|recip", colnames(data_stan_re$data_stan$X_choice))])

# install.packages("corrplot")

color_palette <- colorRampPalette(c("blue", "red"))(100)
relevant_vars_idx <- (rowSums(corr_matrix) > 0) | (colSums(corr_matrix) > 0)
matrix_to_plot <- corr_matrix[relevant_vars_idx, relevant_vars_idx]

png("corr_matrix.png", width = 14, height = 14, units = "in", res = 300)

corrplot(
    corr_matrix,
    method = "color",
    col = color_palette,
    type = "upper",
    addgrid.col = "grey",
    is.corr = TRUE,
    cl.lim = c(0, 60),
    tl.col = "black",
    tl.srt = 90,        # Rotate labels
    tl.cex = 0.9,       # Adjust label font size
    addCoef.col = "black",
    number.cex = 0.01,
    title = CORRELATION_THRESHOLD,
    mar = c(0, 0, 1, 0)
)
dev.off()

