install.packages('here')
install.packages(c("remotes"))
remotes::install_github("stocnet/goldfish@develop") # Run only once
# remotes::install_github("snlab-ch/goldfish.latent") # Run only once
# remotes::install_version("stocnet/goldfish")
remotes::install_github('snlab-ch/goldfish.latent@develop') # Run only once
# remotes::install_github('https://github.com/snlab-ch/goldfish.latent.git')
# remotes::install_local('/home/milky/drosophila-pipelnie/6_dynam_models/goldfish.latent-develop')
library(goldfish.latent)
library(goldfish)
packageVersion("goldfish")
library(here) # For robust file paths
library(glue) # For string formatting
library(cmdstanr)
source(here::here("src/utils.R"))

# library(foreach)
# library(doParallel)
# n_cores <- parallel::detectCores() - 2 
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)

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

all_groups_data_list <- list()
all_groups_info <- data.frame()
group_counter <- 1
all_groups_cox_coeff_list <- list()
all_groups_sample_cox_coeff_list <- list()
all_groups_samples <- list()

##! add foreach for paralelisation
for (treatment_name in treatments_list) {
    folder_path <- file.path(edgelist_path, treatment_name)
    csv_edgelists <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
    # foreach
    for (file_path in csv_edgelists) {
        print(group_name)
        group_name <- tools::file_path_sans_ext(basename(file_path))
        interaction_data <- read.csv(file_path)
        interaction_data <- interaction_data[, c('time', 'sender', 'receiver', 'increment')]
        
        flies <- getActors(interaction_data)
        n_flies <- nrow(flies)
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

        if (treatment_name == "CsCh") {flies$young <- 1}
        else if (treatment_name == "CS_10D") {flies$old <- 1}
        else if (treatment_name == "Cs_5DIZ") {flies$isolated <- 1}

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

        nodesAttr <- link_events(x = nodesAttr, change_event = cov_in_degree, attribute = 'popularity')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_out_degree, attribute = 'activity')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_in_weighted, attribute = 'popularity_weighted')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_out_weighted, attribute = 'activity_weighted')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_positiveinfluence, attribute = 'positive_influence')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_negativeinfluence, attribute = 'negative_influence')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_positiveinf_weighted, attribute = 'positive_inf_weighted')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_negativeinf_weighted, attribute = 'negative_inf_weighted')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_distance_traveled_between_interactions, attribute = 'distance_traveled_between_interactions')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_number_of_flies_in_soc_space, attribute = 'number_of_flies_in_soc_space')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_unique_partners_met_interaction_space, attribute = 'unique_partners_met_interaction_space')
        nodesAttr <- link_events(x = nodesAttr, change_event = cov_unique_partners_met_social_space, attribute = 'unique_partners_met_social_space')

        group_ssm_path <- file.path(soc_space_matrix_path, treatment_name, group_name)
        initial_social_space_matrix <- as.matrix(read.csv(file.path(group_ssm_path, "m_at_t0.csv"), row.names = 1))
        social_space_updates <- read.csv(file.path(group_ssm_path, "updates.csv"))

        if (colnames(social_space_updates)[1] == "X") social_space_updates <- social_space_updates[, -1]
        names(social_space_updates) <- c("time", "sender", "receiver", "replace")
        social_space_updates$sender <- glue("fly{sender+1}",sender=social_space_updates$sender)
        social_space_updates$receiver <- glue("fly{receiver+1}",receiver=social_space_updates$receiver)

        correct_order <- flies$label
        initial_social_space_matrix <- initial_social_space_matrix[correct_order, correct_order]

        socialSpace_network <- make_network(
            matrix = initial_social_space_matrix, 
            nodes = nodesAttr, 
            directed = FALSE
        )
        socialSpace_network <- link_events(
            socialSpace_network, 
            change_event = social_space_updates, 
            nodes = nodesAttr
        )

        dataDynam <- make_data(nodesAttr, interaction_network, dependent, socialSpace_network, flies)

        data_for_this_group_choice <- make_data_re(
            random_effects = list(inertia(interaction_network, weighted = FALSE, window = 288) ~ 1), #! WHAT RE SHOULD WE USE
            fixed_effects = dependent ~
                                indeg(interaction_network, weighted = TRUE, window = 96) +
                                indeg(interaction_network, weighted = TRUE, window = 288) +
                                indeg(interaction_network, weighted = TRUE, window = 864) +
                                indeg(interaction_network, weighted = TRUE, window = 2592) +

                                outdeg(interaction_network, weighted = TRUE, window = 96) +
                                outdeg(interaction_network, weighted = TRUE, window = 288) +
                                outdeg(interaction_network, weighted = TRUE, window = 864) +
                                outdeg(interaction_network, weighted = TRUE, window = 2592) +

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
                                alter(flies$isolated),

            support_constraint = ~ tie(socialSpace_network),
            model = "DyNAM",
            sub_model = "choice",
            data = dataDynam
        )
        
        all_groups_data_list[[group_name]] <- data_for_this_group_choice
        all_groups_info <- rbind(all_groups_info, data.frame(
            group = group_name,
            treatment = treatment_name,
            ixGroup = group_counter
        ))
        group_counter <- group_counter + 1

        correlation_matrix <- cor(data_for_this_group_choice$dataStan$Xchoice)
        high_corr_pairs <- which(abs(correlation_matrix) > 0.8 & correlation_matrix != 1, arr.ind = TRUE)
        if (nrow(high_corr_pairs) > 0) {
            log_file <- paste0("high_correlation_warnings_", group_name, ".txt")
            sink(log_file)
            unique_pairs <- high_corr_pairs[high_corr_pairs[,1] < high_corr_pairs[,2], , drop = FALSE]
            if (nrow(unique_pairs) > 0) {
                print(unique_pairs)
            }
            sink()  # Stop writing to file
        }

        samplesGroups <- sample_data(data_for_this_group_choice, sample_size_set=1)
        all_groups_samples[[group_name]] <- samplesGroups
    }
}

saveRDS(all_groups_data_list, file = "all_groups_data_list.rds")
saveRDS(all_groups_samples, file = "all_groups_samples.rds")



###? -------------------------------------------------------------------------------------------------------------------

        dataCox <- make_cox_model_data(data_for_this_group_choice, "choice")
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


    
        # data_for_this_group_rate <- CreateData(
        #     randomEffects = list(indeg(interaction_network, weighted = TRUE, window = 288) ~ 1),
        #     fixedEffects = dependent ~ 1 +
        #                         indeg(interaction_network, weighted = TRUE, window = 96) +
        #                         # indeg(interaction_network, weighted = TRUE, window = 288) +
        #                         indeg(interaction_network, weighted = TRUE, window = 864) +
        #                         indeg(interaction_network, weighted = TRUE, window = 2592) +
            
        #                         outdeg(interaction_network, weighted = TRUE, window = 96) +
        #                         outdeg(interaction_network, weighted = TRUE, window = 288) +
        #                         outdeg(interaction_network, weighted = TRUE, window = 864) +
        #                         outdeg(interaction_network, weighted = TRUE, window = 2592) +

        #                         indeg(interaction_network, weighted = FALSE, window = 96) +
        #                         indeg(interaction_network, weighted = FALSE, window = 288) +
        #                         indeg(interaction_network, weighted = FALSE, window = 864) +
        #                         indeg(interaction_network, weighted = FALSE, window = 2592) +

        #                         outdeg(interaction_network, weighted = FALSE, window = 96) +
        #                         outdeg(interaction_network, weighted = FALSE, window = 288) +
        #                         outdeg(interaction_network, weighted = FALSE, window = 864) +
        #                         outdeg(interaction_network, weighted = FALSE, window = 2592) +

        #                         nodeTrans(interaction_network, window = 96) +
        #                         nodeTrans(interaction_network, window = 288) +
        #                         nodeTrans(interaction_network, window = 864) +
        #                         nodeTrans(interaction_network, window = 2592) +

        #                         ego(nodesAttr$positive_influence) +
        #                         ego(nodesAttr$negative_influence) +
        #                         ego(nodesAttr$positive_inf_weighted) +
        #                         ego(nodesAttr$negative_inf_weighted) +

        #                         ego(nodesAttr$distance_traveled_between_interactions) + #! difference in distance traveled between interactions 
        #                         ego(nodesAttr$number_of_flies_in_soc_space) + #! number of flies in soc space
        #                         ego(nodesAttr$unique_partners_met_interaction_space) + #! number of uniques fly met between interaction
        #                         ego(nodesAttr$unique_partners_met_social_space) + #! number of uniques fly met between interaction

        #                         ego(flies$young) +
        #                         ego(flies$old) +
        #                         ego(flies$isolated),
        #     supportConstraints = ~ outdeg(socialSpace_network, transformFun = \(x) ifelse(x >= 1, 1, 0)),
        #     model = "DyNAM",
        #     subModel = "rate"
        # )
#     }
# }



# save data_for_this_group with saveRDS .rds
#####
# read the data_for_this_group






library(survival)
library(tidy)

data_stan_re <- GatherGroups(
    groupsInfo = all_groups_info,
    data = all_groups_data_list
)

stanCode <- CreateModelCode(data_stan_re)
modRandomEffects <- cmdstan_model(stanCode)

random_effects_model <- modRandomEffects$sample(
    x = data_stan_re,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 200 #1000,
    iter_sampling = 300 #1500
)

print(summary(random_effects_model))
model_save_path <- file.path(results_path, "full_random_effects_model.rds")
saveRDS(random_effects_model, file = model_save_path)
print(glue("Full random effects model object saved to: {model_save_path}"))


