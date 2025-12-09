remotes::install_github("snlab-ch/goldfish.latent@develop")

check98 <- make_data_re(
  random_effects = NULL,
  fixed_effects = callsDependent ~ 1 + indeg + outdeg,
  sub_model = "rate",
  data = socialEvolutionData,
  support_constraint = NULL,
  progress = TRUE
)

timesInt <- ceiling(check98$data_stan$T_rate / 3)
check98$data_stan$interaction <- rep(1:3, each = timesInt, length.out = check98$data_stan$T_rate)
check98$data_stan$C <- 3

check99 = write_json(
  check98,
  file_name = "../extra/check99.json",
  grain_size = 146
)

stan_code <- make_model_code(check98)

mod_rate <- cmdstan_model("DNRE_Q0_rate_interaction_map_reduce.stan", 
                          cpp_options = list(stan_threads = TRUE))

mod_rate_samples <- mod_rate$sample(
  data = "../extra/check99.json",
  parallel_chains = 1, chains = 1,  iter_warmup = 100, iter_sampling = 100,
  threads_per_chain = 3
)



function (x, file_name, n_chunks = 10, grain_size = 5) 
{
    model <- attr(x, "model")
    sub_model <- attr(x, "sub_model")
    has_sample <- attr(x, "sample")
    scale <- attr(x, "scale")
    data_gathered <- x[!grepl("data_stan", names(x))]
    data_gathered$json_file <- file
    class(data_gathered) <- class(x)
    attr(data_gathered, "model") <- model
    attr(data_gathered, "subModel") <- sub_model
    attr(data_gathered, "sample") <- has_sample
    attr(data_gathered, "scale") <- scale
    attr(data_gathered, "json_file") <- TRUE
    data_stan <- x$data_stan
    n_size <- data_stan[[glue("N_{sub_model}")]]
    t_size <- data_stan[[glue("T_{sub_model}")]]
    p_size <- data_stan[[glue("P_{sub_model}")]]
    q_size <- data_stan[[glue("Q_{sub_model}")]]
    x_name <- glue("X_{sub_model}")
    z_name <- glue("Z_{sub_model}")
    has_z <- q_size > 0
    has_interaction <- !is.null(data_stan[["interaction"]])
    is_integer_x <- rlang::is_integerish(data_stan[[x_name]])
    approx_size_x <- approx_nchar_matrix(n_size, p_size, max_val = max(data_stan[[x_name]]), 
        is_integer = is_integer_x)
    if (has_z) {
        is_integer_z <- rlang::is_integerish(data_stan[[z_name]])
        approx_size_z <- approx_nchar_matrix(n_size, q_size, 
            max_val = max(data_stan[[z_name]]), is_integer = is_integer_z)
    }
    else {
        approx_size_z <- 0
        data_stan[[z_name]] <- NULL
    }
    if (has_interaction) {
        approx_size_interaction <- approx_nchar_vector(t_size, 
            100)
    }
    else {
        approx_size_interaction <- 0
    }
    if (sub_model == "rate") {
        mode(data_stan[["is_dependent"]]) <- "integer"
        approx_size_is_dependent <- approx_nchar_vector(t_size, 
            1)
        is_integer_timespan <- rlang::is_integerish(data_stan[["timespan"]])
        approx_size_timespan <- approx_nchar_vector(t_size, max(data_stan[["timespan"]]), 
            is_integer = is_integer_timespan)
        approx_size_rate <- approx_size_is_dependent + approx_size_timespan
    }
    else {
        approx_size_rate <- 0
    }
    approx_max_size <- max(approx_size_x, approx_size_z)
    aprox_size <- approx_nchar_vector(t_size, n_size)
    total_size <- aprox_size * 3 + approx_size_x + approx_size_z + 
        approx_size_rate + approx_size_interaction + 200
    if (total_size < (2^31 - 1)) {
        data_stan[["grain_size"]] <- grain_size
        cmdstanr::write_stan_json(data = data_stan, file = file_name)
        return(data_gathered)
    }
    if (approx_max_size > (2^31 - 1)) {
        row_partition <- parallel::splitIndices(n_size, n_chunks)
        x_text <- c(glue("\"X_{sub_model}\": ["), vapply(row_partition, 
            create_json_chunk_matrix, character(1), matrix = data_stan[[x_name]]), 
            "],")
        stringr::str_sub(x_text[n_chunks + 1], -1) <- ""
        if (has_z) {
            z_text <- c(glue("\"Z_{sub_model}\": ["), vapply(row_partition, 
                create_json_chunk_matrix, character(1), matrix = data_stan[[z_name]]), 
                "],")
            stringr::str_sub(z_text[n_chunks + 1], -1) <- ""
        }
    }
    else {
        x_text <- c(glue("\"X_{sub_model}\":"), toJSON(data_stan[[x_name]], 
            pretty = TRUE, digits = NA), ",")
        if (has_z) {
            z_text <- c(glue("\"Z_{sub_model}\":"), toJSON(data_stan[[z_name]], 
                pretty = TRUE, digits = NA), ",")
        }
    }
    if (aprox_size > (2^31 - 1)) {
        vector_partition <- parallel::splitIndices(t_size, n_chunks)
        start_text <- c(glue("\"start_{sub_model}\": ["), vapply(vector_partition, 
            create_json_chunk_vector, character(1), vector = data_stan[[glue("start_{sub_model}")]]), 
            "],")
        stringr::str_sub(start_text[n_chunks + 1], -1) <- ""
        end_text <- c(glue("\"end_{sub_model}\": ["), vapply(vector_partition, 
            create_json_chunk_vector, character(1), vector = data_stan[[glue("end_{sub_model}")]]), 
            "],")
        stringr::str_sub(end_text[n_chunks + 1], -1) <- ""
        chose_text <- c(glue("\"chose_{sub_model}\": ["), vapply(vector_partition, 
            create_json_chunk_vector, character(1), vector = data_stan[[glue("chose_{sub_model}")]]), 
            "],")
        stringr::str_sub(chose_text[n_chunks + 1], -1) <- ""
        sender_text <- c("\"sender\": [", vapply(row_partition, 
            create_json_chunk_vector, character(1), vector = data_stan[["sender"]]), 
            "],")
        stringr::str_sub(sender_text[n_chunks + 1], -1) <- ""
        if (has_interaction) {
            interaction_text <- c("\"interaction\": [", vapply(vector_partition, 
                create_json_chunk_vector, character(1), vector = data_stan[["interaction"]]), 
                ",")
            stringr::str_sub(interaction_text[n_chunks + 1], 
                -1) <- ""
        }
        if (sub_model == "rate") {
            timespan_text <- c("\"timespan\": [", vapply(vector_partition, 
                create_json_chunk_vector, character(1), vector = data_stan[["timespan"]]), 
                "],")
            stringr::str_sub(timespan_text[n_chunks + 1], -1) <- ""
            is_dependent_text <- c("\"is_dependent\": [", vapply(vector_partition, 
                create_json_chunk_vector, character(1), vector = data_stan[["is_dependent"]]), 
                "],")
            stringr::str_sub(is_dependent_text[n_chunks + 1], 
                -1) <- ""
        }
    }
    else {
        start_text <- c(glue("\"start_{sub_model}\":"), toJSON(data_stan[[glue("start_{sub_model}")]], 
            pretty = TRUE), ",")
        end_text <- c(glue("\"end_{sub_model}\":"), toJSON(data_stan[[glue("end_{sub_model}")]], 
            pretty = TRUE), ",")
        chose_text <- c(glue("\"chose_{sub_model}\":"), toJSON(data_stan[[glue("chose_{sub_model}")]], 
            pretty = TRUE), ",")
        sender_text <- c("\"sender\":", toJSON(data_stan[["sender"]], 
            pretty = TRUE), ",")
        if (has_interaction) {
            interaction_text <- c("\"interaction\":", toJSON(data_stan[["interaction"]], 
                pretty = TRUE), ",")
        }
        if (sub_model == "rate") {
            timespan_text <- c("\"timespan\":", toJSON(data_stan[["timespan"]], 
                pretty = TRUE), ",")
            is_dependent_text <- c("\"is_dependent\":", toJSON(data_stan[["is_dependent"]], 
                pretty = TRUE), ",")
        }
    }
    keep_dttxt <- c(glue("{data}_{sub_model}", data = c("N", 
        "T", "P", "Q")), c("A", if (has_interaction) "C" else NULL))
    data_text <- transform_json(c(data_stan[keep_dttxt], list(grain_size = grain_size)), 
        first_position = 3, first_replace = "", last_replace = "}")
    conn <- file(file_name, open = "wb")
    writeLines("{", conn, useBytes = TRUE)
    writeLines(start_text, conn)
    writeLines(end_text, conn)
    writeLines(chose_text, conn)
    writeLines(sender_text, conn)
    writeLines(x_text, conn)
    if (has_z) {
        writeLines(z_text, conn)
    }
    if (has_interaction) {
        writeLines(interaction_text, conn)
    }
    if (sub_model == "rate") {
        writeLines(timespan_text, conn)
        writeLines(is_dependent_text, conn)
    }
    writeLines(data_text, conn)
    close(conn)
    data_gathered
}