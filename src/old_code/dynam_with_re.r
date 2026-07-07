#%%
# install.packages("jsonlite")
# install.packages("cmdstanr")
remotes::install_github("snlab-ch/goldfish.latent")
library(goldfish.latent)
library(goldfish)
library(jsonlite)
library(dplyr)
library(stringr)
library(here)  # # helps with path specification NOT REQUIRED
library(glue)  # # helps working with string literals NOT REQUIRED
library(cmdstanr)

source(here::here("src/utils.R"))

treatments_list <- c('CS_10D', 'Cs_5DIZ', 'CsCh')

base_path <- "/home/milky/drosophila-pipelnie/6_dynam_models"
edgelist_path <- here(base_path, "data", "edgelists")
covariances_path <- here(base_path, "data", "covariances")
results_path <- here(base_path, "res")

for (treatment_name in treatments_list){
    folder_path <- paste0(edgelist_paht, treatment_name, "/")
    csv_edgelists <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

    for (file_path in csv_edgelists){
        df <- read.csv(file_path) # csv_edgelists
        interaction_data <- df[, c('time', 'sender', 'receiver', 'increment')]
        flies <- getActors(interaction_data)
        flies$bursty <- rep(0,12)
        flies$influence_pos <- rep(0,12)
        flies$influence_neg <- rep(0,12)
        flies$activity <- rep(0,12)
        flies$popularity <- rep(0,12)

        group_name <- basename(file_path)
        group_name <- gsub("\\.csv$", "", group_name)
        group_path = paste0(covariances_path,'/',treatment_name ,'/', group_name , '/')

        cov_burstines_out_degree <- read_and_clean('burstines_out_degree.csv')
        cov_burstines_out_weighted <- read_and_clean('burstines_out_weighted.csv')
        cov_activity <- read_and_clean('in_degree.csv')
        cov_activity <- read_and_clean('out_degree.csv')
        cov_activity_weighted <- read_and_clean('in_weighted.csv')
        cov_activity_weighted <- read_and_clean('out_weighted.csv')
        cov_negativeinfluence <- read_and_clean('negativeinfluence.csv')
        cov_negativeinf_weighted <- read_and_clean('negativeinf_weighted.csv')
        cov_positiveinfluence <- read_and_clean('positiveinfluence.csv')
        cov_positiveinf_weighted <- read_and_clean('positiveinf_weighted.csv')

        nodesAttr <- defineNodes(flies)
        interaction_network <- defineNetwork(nodes = nodesAttr, directed = TRUE)
        interaction_network <- linkEvents(x = interaction_network, changeEvent = interaction_data, nodes = nodesAttr)
        dependent <- defineDependentEvents(events = interaction_data, nodes = nodesAttr,
                                            defaultNetwork = interaction_network)
        
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_burstines_out_degree = 'burstines_out_degree' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_burstines_out_weighted = 'burstines_out_weighted' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_in_degree = 'popularity' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_out_degree = 'activity' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_in_weighted = 'popularity_weighted' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_out_weighted = 'activity_weighted' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_positiveinfluence = 'positive_influence' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_negativeinfluence = 'negative_influence' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_positiveinf_weighted = 'positive_inf_weighted' )
        nodesAttr <- linkEvents(x = nodesAttr, changeEvent = cov_negativeinf_weighted = 'negative_inf_weighted' )

        # Rate submodel
        rate_model <- estimate(dependent ~ 1
                                + indeg(interaction_network, weighted = TRUE, window = 288)
                                + indeg(interaction_network, weighted = TRUE, window = 864)
                                + indeg(interaction_network, weighted = TRUE, window = 2592)
                                + outdeg(interaction_network, weighted = TRUE, window = 864)
                                + outdeg(interaction_network, weighted = TRUE, window = 288)
                                + outdeg(interaction_network, weighted = TRUE, window = 2592)
                                + ego(nodesAttr$bursty)
                                + ego(nodesAttr$positive_influence)
                                + ego(nodesAttr$negative_influence)
                                + ego(nodesAttr$positive_inf_weighted)
                                + ego(nodesAttr$negative_inf_weighted),
                                model = "DyNAM", subModel = "rate")
        
        result_element <- summary(rate_model)
        json_string <- convert_to_json(result_element)
        save_path <- paste0(results_rate_model_path, group_name, '.json')
        writeLines(json_string, save_path)


        dataGroup <- CreateData(
            randomEffects = list(inertia ~ 1),
            fixedEffects = dependent ~ indeg(interaction_network, weighted = TRUE, window = 288)
                            + indeg(interaction_network, weighted = TRUE, window = 864)
                            + indeg(interaction_network, weighted = TRUE, window = 2592)
                            + outdeg(interaction_network, weighted = TRUE, window = 864)
                            + outdeg(interaction_network, weighted = TRUE, window = 288)
                            + outdeg(interaction_network, weighted = TRUE, window = 2592)
                            + ego(nodesAttr$positive_influence)
                            + ego(nodesAttr$negative_influence)
                            + ego(nodesAttr$positive_inf_weighted)
                            + ego(nodesAttr$negative_inf_weighted)
        )
        # + ego(nodesAttr$bursty)
        str(dataGroup)
        cor(dataGroup$dataStan$Xchoice)

        dataSample <- SampleData(
            dataGroup, fractionChoiceSet = 1/11, methodChoiceSet = "srswor",
            methodEvents = "srswor"
        )
        groupInfo <- data.frame(
            group = sprintf("group%02d", seq_len(5)),
            ixGroup = seq_len(5)
        )

        data <- list()
        for (gr in seq_len(5)) {
            data[[gr]] <- dataSample
        }
        names(data) <- sprintf("group%02d", seq_len(5))

        # gather the data all groups
        data1REStan <- GatherGroups(
            groupsInfo = groupInfo,
            data = data
        ) 

        cmdstanr::write_stan_json(
            data1REStan[!names(data1REStan) %in% "groupEquivalence"],
            file = here("input", glue("dtStanMod{modVersion}Version{verOutput}.json"))
        )



        
        # Choice submodel 
        choice_model <- estimate(dependent ~
                                    inertia
                                + recip(interaction_network, weighted = FALSE, window = 288)
                                + recip(interaction_network, weighted = FALSE, window = 864)
                                + recip(interaction_network, weighted = FALSE, window = 2592)
                                + indeg(interaction_network, weighted = TRUE, window = 288)
                                + indeg(interaction_network, weighted = TRUE, window = 864)
                                + indeg(interaction_network, weighted = TRUE, window = 2592)
                                + outdeg(interaction_network, weighted = TRUE, window = 288)
                                + outdeg(interaction_network, weighted = TRUE, window = 864)
                                + outdeg(interaction_network, weighted = TRUE, window = 2592)
                                + trans(interaction_network, weighted = FALSE, window = 288)
                                + trans(interaction_network, weighted = FALSE, window = 864)
                                + trans(interaction_network, weighted = FALSE, window = 2592)
                                + sim(nodesAttr$activity)
                                + sim(nodesAttr$popularity),
                                model = "DyNAM", subModel = "choice",
                                estimationInit = list(returnIntervalLogL = TRUE, maxIterations = 1000))
        
        result_element <- summary(choice_model)
        json_string <- convert_to_json(result_element)
        save_path <- paste0(results_choice_model_path, group_name, '.json')
        writeLines(json_string, save_path)
        print(group_name)
    }
}
