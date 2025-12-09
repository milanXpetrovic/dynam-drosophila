function (groups_info, data, keep_x = NULL, scale = FALSE) 
{
    if (!is.list(data)) {
        cli::cli_abort(c("{.arg data} must be a list.", x = "You've supplied a {.cls {class(data)}}."))
    }
    if (!is.data.frame(groups_info)) {
        cli::cli_abort(c("{.arg groups_info} must be a data frame.", 
            x = "You've supplied a {.cls {class(groups_info)}}."))
    }
    if (!is.null(keep_x) && !is.character(keep_x)) {
        cli::cli_abort(c("{.arg keep_x} must be a character vector or NULL.", 
            x = "You've supplied a {.cls {class(keep_x)}}."))
    }
    if (length(data) != nrow(groups_info)) {
        cli::cli_abort(c("Length of {.arg data} must match number of rows in {.arg groups_info}.", 
            x = "Length of {.arg data} is {.val {length(data)}},", 
            i = "but {.arg groups_info} has {.val {nrow(groups_info)}} rows."))
    }
    model <- attr(data[[1]], "model")
    sub_model <- attr(data[[1]], "sub_model")
    has_sample <- attr(data[[1]], "sample")
    first_group_data <- data[[1]]$data_stan
    x_name <- paste0("X_", sub_model)
    z_name <- paste0("Z_", sub_model)
    n_name <- paste0("N_", sub_model)
    t_name <- paste0("T_", sub_model)
    col_names_x <- colnames(first_group_data[[x_name]])
    n_size <- sum(sapply(data, function(x) x$data_stan[[n_name]]))
    t_size <- sum(sapply(data, function(x) x$data_stan[[t_name]]))
    p_size <- length(col_names_x)
    if (!is.null(keep_x)) {
        if (!all(keep_x %in% col_names_x)) 
            cli::cli_abort(c("Some variables in {.arg keep_x} are not in the design matrix.", 
                x = "You've supplied {setdiff(keep_x, col_names_x)}."))
        p_size <- length(keep_x)
        col_names_x <- keep_x
    }
    q_size <- max(sapply(data, function(x) x$data_stan[[paste0("Q_", 
        sub_model)]]))
    has_z <- !is.null(first_group_data[[z_name]])
    x_matrix <- matrix(0, nrow = n_size, ncol = p_size, dimnames = list(NULL, 
        col_names_x))
    if (has_z) {
        z_matrix <- matrix(0, nrow = n_size, ncol = q_size, dimnames = list(NULL, 
            colnames(first_group_data[[z_name]])))
    }
    n_groups <- nrow(groups_info)
    start_agg <- integer(t_size)
    end_agg <- integer(t_size)
    chose_agg <- integer(t_size)
    sender_agg <- integer(n_size)
    start_group_agg <- integer(n_groups)
    if (sub_model == "rate") {
        timespan_agg <- integer(t_size)
        is_dependent_agg <- integer(t_size)
    }
    n_offset <- 0
    t_offset <- 0
    for (gr in seq_len(n_groups)) {
        group_data <- data[[gr]]$data_stan
        group_label <- groups_info$ixGroup[gr]
        current_col_names <- colnames(group_data[[x_name]])
        if (!is.null(keep_x)) {
            missing_names <- setdiff(keep_x, current_col_names)
            if (length(missing_names) > 0) {
                cli::cli_abort(c("Group {group_label} has missing variables defined in", 
                  "{.var keep_x}.", x = "Missing: {missing_names}"))
            }
        }
        else {
            if (!identical(col_names_x, current_col_names)) {
                not_first <- setdiff(col_names_x, current_col_names)
                not_current <- setdiff(current_col_names, col_names_x)
                cli::cli_abort(c("Variables do not match among groups:", 
                  i = "In first group but not in current: {not_first}", 
                  i = "In current but not in first: {not_current}"))
            }
        }
        x_aux <- group_data[[x_name]][, col_names_x, drop = FALSE]
        n_indices <- (n_offset + 1):(n_offset + group_data[[n_name]])
        t_indices <- (t_offset + 1):(t_offset + group_data[[t_name]])
        x_matrix[n_indices, ] <- x_aux
        if (has_z) 
            z_matrix[n_indices, ] <- group_data[[z_name]]
        start_agg[t_indices] <- group_data[[paste0("start_", 
            sub_model)]] + n_offset
        end_agg[t_indices] <- group_data[[paste0("end_", sub_model)]] + 
            n_offset
        chose_agg[t_indices] <- group_data[[paste0("chose_", 
            sub_model)]] + n_offset
        if (sub_model == "rate") {
            timespan_agg[t_indices] <- group_data$timespan
            is_dependent_agg[t_indices] <- group_data$is_dependent
        }
        sender_agg[n_indices] <- groups_info$ixGroup[gr]
        start_group_agg[gr] <- t_indices[1]
        n_offset <- n_offset + group_data[[n_name]]
        t_offset <- t_offset + group_data[[t_name]]
    }
    data_stan <- list(T = t_size, N = n_size, P = p_size, Q = q_size, 
        start = start_agg, end = end_agg, X = x_matrix, Z = if (has_z) z_matrix else NULL, 
        chose = chose_agg)
    names(data_stan) <- glue("{name}_{sub_model}", name = names(data_stan))
    data_stan[["A"]] <- n_groups
    data_stan[["sender"]] <- sender_agg
    data_stan[["start_group"]] <- start_group_agg
    if (sub_model == "rate") {
        data_stan[["timespan"]] <- timespan_agg
        data_stan[["is_dependent"]] <- is_dependent_agg
        mean_actors <- mean(end_agg - start_agg + 1)
        total_time <- sum(timespan_agg)
        data_stan[["offset_int_rate"]] <- log(t_size/(total_time * 
            mean_actors))
    }
    data_gathered <- list(data_stan = data_stan, namesEffects = data[[1]]$namesEffects, 
        effectsDescription = data[[1]]$effectsDescription, groupInfo = groups_info)
    class(data_gathered) <- class(data[[1]])
    attr(data_gathered, "model") <- model
    attr(data_gathered, "subModel") <- sub_model
    attr(data_gathered, "sample") <- has_sample
    return(data_gathered)
}