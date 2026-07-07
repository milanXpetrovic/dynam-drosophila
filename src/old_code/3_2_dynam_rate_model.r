remotes::install_github('snlab-ch/goldfish.latent@develop')
install.packages("lme4")

library(goldfish.latent)
library(goldfish)
library(here)
library(glue)
library(cmdstanr)
library(foreach)
library(doParallel)
library(parallel)
library(iterators)

library(lme4)
library(broom)

source(here::here("src/utils.R"))
treatments_list <- c('CS_10D', 'Cs_5DIZ', 'CsCh')

base_path <- "/home/milky/drosophila-pipelnie/6_dynam_models"
edgelist_path <- file.path(base_path, "data", "edgelists")
covariances_path <- file.path(base_path, "data", "covariances")
soc_space_matrix_path <- file.path(base_path, "data", "soc_space_matrix")
results_path <- file.path(base_path, "res")
results_rate_model_path <- file.path(results_path, "rate_model_individual")

dir.create(results_rate_model_path, showWarnings = FALSE, recursive = TRUE)

log_file_path <- file.path(results_path, "parallel_execution_log.txt")

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

rate_model_formula <- dependent ~ 1 +
                        # indeg(interaction_network, weighted = TRUE, window = 288) +
                        # indeg(interaction_network, weighted = TRUE, window = 864) +
                        # indeg(interaction_network, weighted = TRUE, window = 2592) +
    
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

                        # node_trans(interaction_network, window = 96) +
                        # node_trans(interaction_network, window = 288) +
                        # node_trans(interaction_network, window = 864) +
                        # node_trans(interaction_network, window = 2592) +

                        ego(nodesAttr$positive_influence) +
                        ego(nodesAttr$negative_influence) +
                        ego(nodesAttr$positive_inf_weighted) +
                        ego(nodesAttr$negative_inf_weighted) +

                        ego(nodesAttr$distance_traveled_between_interactions) + #! difference in distance traveled between interactions 
                        ego(nodesAttr$number_of_flies_in_soc_space) + #! number of flies in soc space
                        ego(nodesAttr$unique_partners_met_interaction_space) + #! number of uniques fly met between interaction
                        ego(nodesAttr$unique_partners_met_social_space) + #! number of uniques fly met between interaction

                        ego(flies$young) +
                        ego(flies$old) +
                        ego(flies$isolated)

files_to_process <- do.call(rbind, all_files_list)
# files_to_process <- head(files_to_process, 4) ##! TAKE CARE 

n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)
parallel_results <- foreach(
    i = 1:nrow(files_to_process),
    .packages = c("goldfish", "goldfish.latent", "glue", "tools"),
    .errorhandling = "pass"
) %dopar% {
    worker_pid <- Sys.getpid()
    start_time <- Sys.time()
    start_message <- glue::glue(
        "{format(start_time, '%Y-%m-%d %H:%M:%S')} [PID: {worker_pid}] START - Processing group: {files_to_process$group_name[i]}"
    )
    cat(paste0(start_message, "\n"), file = log_file_path, append = TRUE)
    treatment_name <- files_to_process$treatment_name[i]
    group_name <- files_to_process$group_name[i]
    file_path <- files_to_process$file_path[i]
    interaction_data <- read.csv(file_path)
    interaction_data <- interaction_data[, c('time', 'sender', 'receiver', 'increment')] 
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
    
    nodesAttr <- make_nodes(flies)
    interaction_network <- make_network(nodes = nodesAttr, directed = TRUE)
    interaction_network <- link_events(x = interaction_network, change_event = interaction_data, nodes = nodesAttr)
    dependent <- make_dependent_events(events = interaction_data, nodes = nodesAttr, default_network = interaction_network)

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
    dataDynam <- make_data(nodesAttr, interaction_network, dependent, flies)
    
    print(glue("Processing group: {group_name}"))
    data_for_this_group_rate <- make_data_re(
        random_effects = NULL,
        fixed_effects = rate_model_formula,
        # support_constraint = ~ tie(socialSpace_network),
        model = "DyNAM",
        sub_model = "rate",
        data = dataDynam,
        control_preprocessing = set_preprocessing_opt(start_time = 0) #! , end_time = 100
    )
    end_time <- Sys.time()
    duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    end_message <- glue::glue(
        "{format(end_time, '%Y-%m-%d %H:%M:%S')} [PID: {worker_pid}] END   - Finished group: {group_name}. Duration: {duration}s"
    )
    cat(paste0(end_message, "\n"), file = log_file_path, append = TRUE)
    list(
        group_name = group_name,
        treatment_name = treatment_name,
        group_index = i,
        data_object = data_for_this_group_rate,
        samples_object = NULL
    )
}
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
        print(paste("An error occurred for a group:", result))
    }
}

all_groups_info <- do.call(rbind, all_groups_info_list)
saveRDS(all_groups_data_list, file = "all_groups_data_list_rate.rds")
saveRDS(all_groups_samples, file = "all_groups_samples_rate.rds")
saveRDS(all_groups_info, file = "all_groups_info_rate.rds")

########################################################################################################################

all_groups_data_list <- readRDS("all_groups_data_list_rate.rds")
all_groups_info <- readRDS("all_groups_info_rate.rds")

CORRELATION_THRESHOLD <- 0.9
all_colnames <- unique(unlist(lapply(all_groups_data_list, function(group) {
    colnames(group$data_stan$X_rate[, -1, drop = FALSE])})))
all_colnames <- sort(all_colnames)
final_count_matrix <- matrix(0,
    nrow = length(all_colnames),
    ncol = length(all_colnames),
    dimnames = list(all_colnames, all_colnames)
)
for (group_name in names(all_groups_data_list)) {
    fixed_effects_matrix <- all_groups_data_list[[group_name]]$data_stan$X_rate[, -1, drop = FALSE]
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
# install.packages("corrplot")
library(corrplot)
library(ggplot2)
library(reshape2)
color_palette <- colorRampPalette(c("white", "red"))(100)
relevant_vars_idx <- (rowSums(final_count_matrix) > 0) | (colSums(final_count_matrix) > 0)
matrix_to_plot <- final_count_matrix[relevant_vars_idx, relevant_vars_idx]
png("high_correlation_plot_rate09.png", width = 14, height = 14, units = "in", res = 300)
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
########################################################################################################################

all_groups_data_list <- head(all_groups_data_list)
all_groups_info <- head(all_groups_info)
all_groups_data_list <- readRDS("all_groups_data_list_rate.rds")
all_groups_info <- readRDS("all_groups_info_rate.rds")

data_stan_re <- gather_groups(
    groups_info = all_groups_info,
    data = all_groups_data_list
)

dataCox <-make_cox_model_data(data_stan_re, "rate")
names(dataCox) <- stringr::str_replace_all(names(dataCox), ":", "__")
whichDistance <- grep("ego_nodesAttrOfdistance_traveled_between_interactions", names(dataCox))
dataCox[, whichDistance] <- log1p(dataCox[, whichDistance]) 
termsAdd <- names(dataCox)[seq.int(2, which(names(dataCox) == "event") - 1)]
termsAdd <- names(dataCox)[c(1:81)]
# glm model
mod00 <- glm(
  as.formula(paste0("selected ~ - 1 + offset(log_timespan) + ",
    paste(termsAdd, collapse = " + ")
  )),
  family = poisson(link = "log"),
  data = dataCox[dataCox$log_timespan > -Inf, ]
)

termsInteract <- names(dataCox)[seq.int(2, which(names(dataCox) == "ego_fliesOfyoung") - 1)]
# termsInteract <- names(dataCox)[c(8, 12)]
termsExpCond <- grep("ego_fliesOf", names(dataCox), value = TRUE)
interaction_formula <- as.formula(paste0(
    "selected ~ -1 + offset(log_timespan) +", 
    paste(termsExpCond, collapse = " + "), "+ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ")"
))
mod01 <- glm(
    interaction_formula,
    family = poisson(link = "log"),
    data = dataCox[dataCox$log_timespan > -Inf, ]
    )
coefTable01 <- tidy(mod01, conf.int = T )
write.csv(coefTable01, file = "tmp/res_rate_01.csv")

termsInteract <- names(dataCox)[seq.int(2, which(names(dataCox) == "ego_fliesOfyoung") - 1)]
excludeTerms <- grepl("^node.trans|W.\\d+$", termsInteract)
termsInteract<- termsInteract[!excludeTerms]
# termsInteract <- names(dataCox)[c(8, 12)]
termsExpCond <- grep("ego_fliesOf", names(dataCox), value = TRUE)
interaction_formula <- as.formula(paste0(
    "selected ~ -1 + offset(log_timespan) +", 
    paste(termsExpCond, collapse = " + "), "+ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ")"
))
mod02 <- glm(
    interaction_formula,
    family = poisson(link = "log"),
    data = dataCox[dataCox$log_timespan > -Inf, ]
    )
coefTable02 <- tidy(mod02, conf.int = F )
write.csv(coefTable02, file = "tmp/res_rate_02.csv")

interaction_formula <- as.formula(paste0(
    "selected ~ -1 + (1|sender) + offset(log_timespan) +", 
    paste(termsExpCond, collapse = " + "), "+ (",
    paste(termsExpCond, collapse = " + "), ") : (",
    paste(termsInteract, collapse = " + "),
    ")"
))
mod03 <- glmer(
    interaction_formula,
    family = poisson(link = "log"),
    data = dataCox[dataCox$log_timespan > -Inf, ]
    )
coefTable03 <- tidy(mod03, conf.int = FALSE)
write.csv(coefTable03, file = "res02.csv")

########################################################################################################################
# remotes::install_github("snlab-ch/goldfish.latent@develop")
remotes::install_github("auzaheta/goldfish.latent@develop", build_vignettes = FALSE)
# install.packages("prettyunits")
library(goldfish.latent)
library(goldfish)
library(here)
library(glue)
library(cmdstanr)
library(foreach)
library(doParallel)
library(parallel)
library(iterators)
source(here::here("src/utils.R"))
# setwd('/home/milky/drosophila-pipelnie/6_dynam_models/src')
library(qs)
library(psych)
# qsave(all_groups_data_list, "all_groups_data_list_rate.qs")
all_groups_data_list <- qread("/home/milky/droso-pipe/6_dynam_models/all_groups_data_list_rate.qs")
all_groups_info <- readRDS("/home/milky/droso-pipe/6_dynam_models/all_groups_info_rate.rds")

data_stan_re <- gather_groups(
    groups_info = all_groups_info,
    data = all_groups_data_list
)

data_stan_re$data_stan$grain_size <- 1621
data_stan_re$data_stan$Z_rate <- numeric(0)

old_X_rate <- data_stan_re$data_stan$X_rate
var_distance <- grep("distance", colnames(old_X_rate), value = T)
old_X_rate[, var_distance] <- log1p(old_X_rate[, var_distance])

data_stan_re$data_stan$X_rate_raw <- old_X_rate[, -1]

data_stan_re$data_stan$X_rate <- NULL
data_stan_re$data_stan$P_rate <- ncol(old_X_rate) - 1
data_stan_re$data_stan$send_int <- rep(1:3, each = 20)

data_stan_re$data_stan$interaction <- c(
  rep(1, data_stan_re$data_stan$start_group[20]),
  rep(2, data_stan_re$data_stan$start_group[40] - data_stan_re$data_stan$start_group[20]),
  rep(3, data_stan_re$data_stan$T_rate - data_stan_re$data_stan$start_group[40])
)
data_stan_re$data_stan$C <- 3
data_stan_re$data_stan$sub_model <- NULL
attr(data_stan_re, "sub_model") <- "rate"

# check99 = goldfish.latent::write_json(
#   data_stan = data_stan_re,
#   file_name = "check99.json",
#   n_chunks = 100,
#   grain_size = 10
# )

data_stan_re$data_stan$Z_rate = NULL

json_file <- write_json(
    data_stan = data_stan_re,
    file_name = "rate.json",
    n_chunks = 100,
    grain_size = 10
)

mod_rate <- cmdstan_model("/home/milky/droso-pipe/6_dynam_models/DNRE_Q0_rate_int_centered.stan", cpp_options = list(stan_threads = TRUE))

# mod_rate_samples <- mod_rate$sample(
#   data = data_stan_re$data_stan,
#   parallel_chains = 4, chains = 4,  iter_warmup = 500, iter_sampling = 500,
#   threads_per_chain = 10
# ) 

# data_stan_export <- list(X_rate = data_stan_re$data_stan$X_rate)
# cmdstanr::write_stan_json(data_stan_re$data_stan, file = "tmp/data_rate.json")
# cmdstanr::write_stan_json(data_stan_export, file = "tmp/data_rate.json")

# data_stan_re$data_stan 

random_effects_model <- mod_rate$sample(
    data = json_file, #
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 11,
    iter_warmup = 500,
    iter_sampling = 500,
    output_dir = "/home/milky/droso-pipe/6_dynam_models/src/samples"
)   

random_effects_model$save_object(file = "tmp/full_random_effects_model_rateFINAL.rds")
arrayDraws <- random_effects_model$draws()
saveRDS(arrayDraws, file = "tmp/arrayDraws_rateFINAL.rds")

draws_beta <- random_effects_model_wsw$draws("beta_rate", format = "draws_matrix")
colnames(draws_beta) <- colnames(data_stan_wsw$X_rate)
beta_intervals_original <- mcmc_intervals_data(draws_beta)
write.csv(beta_intervals_original, file = "tmp/beta_intervals_original_rateFINAL.csv")

gamma_intervals <- mcmc_intervals_data(draws_gamma)
write.csv(gamma_intervals, file = "tmp/gamma_intervals_rateFINAL.csv")


###################

get_stats <- function(x) {
  c(
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    min    = min(x, na.rm = TRUE),
    `5%`   = quantile(x, 0.05, na.rm = TRUE),
    `10%`  = quantile(x, 0.10, na.rm = TRUE),
    `25%`  = quantile(x, 0.25, na.rm = TRUE),
    `50%`  = quantile(x, 0.50, na.rm = TRUE),
    `75%`  = quantile(x, 0.75, na.rm = TRUE),
    `90%`  = quantile(x, 0.90, na.rm = TRUE),
    `95%`  = quantile(x, 0.95, na.rm = TRUE),
    max    = max(x, na.rm = TRUE)
  )
}

desc_table <- apply(new_X_rate, 2, get_stats)
desc_table <- t(desc_table)
write.csv(desc_table, "new_X_rate_stats.csv", row.names = TRUE)

install.packages("openxlsx")
library(openxlsx)
write.xlsx(desc_table, file = "new_X_rate_stats.csv.xlsx", rowNames = TRUE)
