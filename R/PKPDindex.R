#' PKPDindex: Optimal PK/PD Index Finder
#'
#' This function fits various Emax models to a given dataset, allowing for flexibility in model selection, initial parameter estimates, and plotting options.
#' @param dataset A data frame containing the independent (x) and dependent (y) variables.
#' @param x_columns A character vector specifying the x-axis variables (PK/PD indices).
#'   If NULL (default), the function attempts to detect appropriate columns from the dataset,
#'   specifically `"auc_mic"`, `"cmax_mic"`, and `"t_mic"`. If these are not found, the user must specify the names manually.
#'   - `"auc_mic"`: area under the concentration-time curve divided by the MIC.
#'   - `"cmax_mic"`: peak drug concentration divided by the MIC.
#'   - `"t_mic"`: time above MIC (duration the drug concentration exceeds MIC).
#'   Users should calculate these indices based on their PK data before using this function.
#' @param y_column A character string specifying the response variable.
#'   Default name is `"response"`.The response should be the log10-transformed change in CFU/ml (Delta log10 CFU/ml).
#'   Users can either provide a column with pre-calculated log10 CFU/ml changes,
#'   or provide raw CFU/ml counts at the initial (CFU_init) and 24-hour timepoint (CFU_24)
#'   , and the function will automatically calculate the log10 change in CFU/ml (Delta log10 CFU/ml).
#' @param E0_fix Fixed E0 (baseline effect) value.
#' @param Emax_fix Fixed Emax (maximum effect) value.
#' @param EI50_init Optional numeric vector specifying initial EI50 values for each x_column. Defaults to NULL, and values are estimated automatically.
#' @param maxiter Maximum number of iterations - Specifies the maximum number of iterations allowed for the nonlinear least squares (NLS) fitting process. Higher values may help convergence for complex models. Default maxiter = 500.
#' @param tol Tolerance level - Defines the tolerance for convergence in the NLS algorithm. Lower values indicate stricter convergence criteria. Default tol = 1e-5.
#' @param minFactor Minimum step factor - Determines the smallest step size used in parameter updates during the NLS fitting process, controlling the precision of optimisation. Default minFactor = 1e-7.
#' @param select_mod Optional named list specifying preferred models for each x_column.
#' @param plot_results Logical; if TRUE, the function generates model fit plots.
#' @param srow Single row plotting - Logical (TRUE or FALSE). If TRUE, plots all best model fits in a single row for visual comparison.
#' @param xlim A numeric vector of length 2 specifying x-axis limits.
#' @param ylim A numeric vector of length 2 specifying y-axis limits.
#' @param point_color Optional character string specifying the point colour in plots.
#' @param line_color Optional character string specifying the line colour in plots.
#' @param x_label Optional named list specifying custom x-axis labels.
#' @param y_label Optional character string specifying a custom y-axis label.
#' @param plot_title Optional character string specifying a custom plot title.
#' @param log_scale_x Optional named list specifying whether to apply log10 scaling to x-axis for each x_column.
#' @param title_cex Size of the plot title text. Default title_cex = 1.2.
#' @param label_cex Size of the axis title. Default label_cex = 1.0.
#' @param axis_cex Size of the axis labels. Default axis_cex = 1.0.
#' @param detail_cex Size of the model detail text on the plot. Default detail_cex = 1.0.
#'
#' @return A list containing:
#'   - **All_Model_Results**: A data frame with results from all fitted models.
#'   - **Best_Models**: A data frame with the best model (lowest AIC) for each PK/PD index.
#'   - **Plots**: A list of recorded plots (if `plot_results = TRUE`).
#' @export
#'
#' @examples
#' # Basic usage with default settings
#'  output <- PKPDindex(
#'   dataset = PKPDindex_data,
#'   E0_fix = 1.5,
#'   Emax_fix = 4.8
#' )
#' # Custom x and y columns and initial data
#'  output <- PKPDindex(
#'   dataset = PKPDindex_data,
#'   E0_fix = 1.5,
#'   Emax_fix = 4.8,
#'   x_columns = c("auc_mic","cmax_mic","t_mic"),
#'   y_column = "response",
#'   EI50_init = c(1,1,1)
#' )
#'
#' # Generate and custom plots
#'  output <- PKPDindex(
#'   dataset = PKPDindex_data,
#'   E0_fix = 1.5,
#'   Emax_fix = 4.8,
#'   plot_results = TRUE,
#'   srow=TRUE,
#'   xlim = c(0, 50),
#'   ylim = c(-2, 10),
#'   point_color = "green",
#'   line_color = "purple",
#'   select_mod = list(auc_mic = "m5", t_mic = "m1"),
#'   x_label = list(auc_mic = "AUC/MIC", cmax_mic = "Cmax/MIC", t_mic = "Time>MIC"),
#'   y_label = "Log10 Change in CFU",
#'   plot_title = "Model Fitting Results",
#'   log_scale_x = list(auc_mic = TRUE, cmax_mic = TRUE,t_mic=FALSE),
#'   title_cex = 2,
#'   label_cex = 1.5,
#'   axis_cex = 1.4,
#'   detail_cex = 1.3
#' )
#'
#' #' # To view the best models:
#' output$Best_Models
#'
#' # To view all model results:
#' output$All_Model_Results
#'
#' # To access a specific plot:
#' output$Plots[["cmax_mic"]]
#' @details
#' The function fits different variations of the Emax model to describe the relationship between PK/PD indices and response.
#' The available models (m1 to m8) are defined as follows:
#'
#' - **m1**: Fixed E0 and Emax, no Hill coefficient.
#' - **m2**: Fixed E0 and Emax, with Hill coefficient (gam).
#' - **m3**: Fixed E0, estimated Emax, no Hill coefficient.
#' - **m4**: Fixed E0, estimated Emax, with Hill coefficient.
#' - **m5**: Estimated E0, fixed Emax, no Hill coefficient.
#' - **m6**: Estimated E0, fixed Emax, with Hill coefficient.
#' - **m7**: Estimated E0 and Emax, no Hill coefficient.
#' - **m8**: Fully estimated model (E0, Emax, EI50, and gam).
#'
#' Users can select specific models using the `select_mod` argument.


#' @importFrom stats AIC as.formula coef nls nls.control residuals
#' @importFrom graphics axis lines par text
#' @importFrom grDevices recordPlot


PKPDindex <- function(dataset, x_columns = NULL, y_column = "response", E0_fix, Emax_fix, EI50_init = NULL ,maxiter = 500, tol = 1e-5, minFactor = 1e-7,select_mod = NULL, plot_results = FALSE ,srow=FALSE, xlim = NULL, ylim = NULL, point_color = NULL,line_color = NULL,  x_label = NULL, y_label = NULL, plot_title = NULL, log_scale_x = NULL,title_cex = 1.2,label_cex = 1.0,axis_cex = 1.0,detail_cex = 1.0) {

  # Convert dataset column names to lowercase for easy matching
  colnames(dataset) <- tolower(colnames(dataset))
  # Handle response or raw CFU data
  if (!(y_column %in% colnames(dataset))) {
    if (all(c("cfu_init", "cfu_24") %in% colnames(dataset))) {
      message("No '", y_column, "' column found. Calculating response from CFU_init and CFU_24.")
      dataset[[y_column]] <- log10(dataset$cfu_24) - log10(dataset$cfu_init)
    } else {
      stop("No '", y_column, "' column found and CFU_init/CFU_24 not available.")
    }
  }

  # Handle default or user-defined y_column
  y_column <- tolower(y_column)

  if (!y_column %in% colnames(dataset)) {
    stop(paste0("The specified y_column '", y_column, "' does not exist in the dataset."))
  }

  # Default PK/PD index candidates
  default_x_columns <- c("auc_mic", "cmax_mic", "t_mic")

  # Auto-detect x_columns if not specified
  if (is.null(x_columns)) {
    detected_columns <- intersect(default_x_columns, colnames(dataset))

    if (length(detected_columns) == 0) {  # Allow any number instead of requiring three
      stop("No valid default columns found. Please specify x_columns.")
    }
    x_columns <- detected_columns
  }

  # Ensure selected columns are numeric
  if (!all(sapply(dataset[, c(x_columns, y_column)], is.numeric))) {
    stop("All selected columns in the dataset must be numeric.")
  }

  # Compute initial EI50 for each PK/PD index or use user-specified values
  if (!is.null(EI50_init)) {
    if (length(EI50_init) != length(x_columns)) {
      stop("The length of EI50_init must match the number of x_columns.")
    }
    EI50_init_values <- EI50_init
  } else {
    EI50_init_values <- sapply(x_columns, function(col) {
      y50_target <- E0_fix - Emax_fix / 2
      closest_point <- which.min(abs(dataset[[y_column]] - y50_target))
      dataset[[col]][closest_point]
    })
  }

  gam_init <- 1

  # Define model formulas dynamically
  model_formulas <- list(
    m1 = function(x) paste0(y_column, " ~ E0_fix - Emax_fix * ", x, " / (EI50 + ", x, ")"),
    m2 = function(x) paste0(y_column, " ~ E0_fix - Emax_fix * ", x, "^gam / (EI50^gam + ", x, "^gam)"),
    m3 = function(x) paste0(y_column, " ~ E0_fix - Emax * ", x, " / (EI50 + ", x, ")"),
    m4 = function(x) paste0(y_column, " ~ E0_fix - Emax * ", x, "^gam / (EI50^gam + ", x, "^gam)"),
    m5 = function(x) paste0(y_column, " ~ E0 - Emax_fix * ", x, " / (EI50 + ", x, ")"),
    m6 = function(x) paste0(y_column, " ~ E0 - Emax_fix * ", x, "^gam / (EI50^gam + ", x, "^gam)"),
    m7 = function(x) paste0(y_column, " ~ E0 - Emax * ", x, " / (EI50 + ", x, ")"),
    m8 = function(x) paste0(y_column, " ~ E0 - Emax * ", x, "^gam / (EI50^gam + ", x, "^gam)")
  )

  all_results <- list()
  plots <- list()
  best_models <- list()

  for (i in seq_along(x_columns)) {
    x_column <- x_columns[i]
    EI50_inital <- EI50_init_values[i]

    model_results <- data.frame(Model = character(0), AIC = numeric(0), R_squared = numeric(0),
                                E0 = numeric(0), Emax = numeric(0), gam = character(0), EI50 = numeric(0), Success = character(0), stringsAsFactors = FALSE)

    # Check if the user specified models for this index
    selected_models <- if (!is.null(select_mod) && x_column %in% names(select_mod)) {
      select_mod[[x_column]]
    } else {
      names(model_formulas)
    }

    for (model_name in selected_models) {
      fit_result <- try({
        model_formula <- as.formula(model_formulas[[model_name]](x_column))

        start_list <- switch(
          model_name,
          m1 = list(EI50 = EI50_inital),
          m2 = list(EI50 = EI50_inital, gam = gam_init),
          m3 = list(Emax = Emax_fix, EI50 = EI50_inital),
          m4 = list(Emax = Emax_fix, EI50 = EI50_inital, gam = gam_init),
          m5 = list(E0 = E0_fix, EI50 = EI50_inital),
          m6 = list(E0 = E0_fix, EI50 = EI50_inital, gam = gam_init),
          m7 = list(E0 = E0_fix, Emax = Emax_fix, EI50 = EI50_inital),
          m8 = list(E0 = E0_fix, Emax = Emax_fix, EI50 = EI50_inital, gam = gam_init)
        )

        fit <- nls(model_formula,
                   start = start_list,
                   data = dataset,
                   control = nls.control(maxiter = maxiter, tol = tol, minFactor = minFactor))

        r_squared <- 1 - sum(residuals(fit)^2) / sum((dataset[[y_column]] - mean(dataset[[y_column]]))^2)
        aic_value <- AIC(fit)

        fitted_coefficients <- coef(fit)
        parameter_values <- c(
          E0 = ifelse("E0" %in% names(fitted_coefficients), fitted_coefficients["E0"], E0_fix),
          Emax = ifelse("Emax" %in% names(fitted_coefficients), fitted_coefficients["Emax"], Emax_fix),
          gam = ifelse("gam" %in% names(fitted_coefficients), fitted_coefficients["gam"], NA),
          EI50 = ifelse("EI50" %in% names(fitted_coefficients), fitted_coefficients["EI50"], EI50_inital)
        )

        model_results <- rbind(model_results, data.frame(
          Model = model_name,
          AIC = round(aic_value,2),
          R_squared = paste0(round(r_squared * 100, 2), "%"),
          E0 = round(parameter_values["E0"], 2),
          Emax = round(parameter_values["Emax"], 2),
          gam = ifelse(is.na(parameter_values["gam"]), 1, round(parameter_values["gam"], 2)),
          EI50 = round(parameter_values["EI50"], 2),
          Success = "Success"
        ))

      }, silent = TRUE)

      if (inherits(fit_result, "try-error")) {
        model_results <- rbind(model_results, data.frame(
          Model = model_name,
          AIC = NA,
          R_squared = NA,
          E0 = NA,
          Emax = NA,
          gam = NA,
          EI50 = NA,
          Success = "Failure"
        ))
      }
    }

    # Find the model with the lowest AIC
    best_model <- model_results[which.min(model_results$AIC), ]
    best_models[[x_column]] <- best_model

    # Store all model results including success/failure status
    all_results <- c(all_results, list(data.frame(
      PKPD_Index = x_column,
      Model = model_results$Model,
      AIC = model_results$AIC,
      R_squared = model_results$R_squared,
      E0 = model_results$E0,
      Emax = model_results$Emax,
      gam = model_results$gam,
      EI50 = model_results$EI50,
      Success = model_results$Success
    )))
  }

  # Plotting best models
  # Plotting best models
  if (plot_results) {
    oldpar <- par(no.readonly = TRUE)  # Save user's current graphics
    on.exit(par(oldpar))               # Restore them when function exits
    if (srow) {
      par(mfrow = c(1, length(best_models)), mar = c(4, 4, 2, 1))  # One row for all plots
    } else {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))  # Default single plot
    }

    # Define default colors
    dot_color <- "blue"
    fit_color <- "red"

    # Allow user to specify custom colors for dots and lines
    if (!is.null(point_color)) {
      dot_color <- point_color
    }
    if (!is.null(line_color)) {
      fit_color <- line_color
    }

    # Allow user to specify custom labels for y-axis and plot title
    c_y_label <- ifelse(is.null(y_label), "Response", y_label)
    c_plot_title <- ifelse(is.null(plot_title), "Best Model Fit", plot_title)

    # Loop through best_models to apply log scale for each plot
    for (index in names(best_models)) {
      best_model <- best_models[[index]]

      if (is.na(best_model$AIC)) {
        next
      }

      best_model_name <- best_model$Model
      model_formula <- as.formula(model_formulas[[best_model_name]](index))

      start_list <- switch(
        best_model_name,
        m1 = list(EI50 = as.numeric(best_model$EI50)),
        m2 = list(EI50 = as.numeric(best_model$EI50), gam = as.numeric(best_model$gam)),
        m3 = list(Emax = as.numeric(best_model$Emax), EI50 = as.numeric(best_model$EI50)),
        m4 = list(Emax = as.numeric(best_model$Emax), EI50 = as.numeric(best_model$EI50), gam = as.numeric(best_model$gam)),
        m5 = list(E0 = as.numeric(best_model$E0), EI50 = as.numeric(best_model$EI50)),
        m6 = list(E0 = as.numeric(best_model$E0), EI50 = as.numeric(best_model$EI50), gam = as.numeric(best_model$gam)),
        m7 = list(E0 = as.numeric(best_model$E0), Emax = as.numeric(best_model$Emax), EI50 = as.numeric(best_model$EI50)),
        m8 = list(E0 = as.numeric(best_model$E0), Emax = as.numeric(best_model$Emax), EI50 = as.numeric(best_model$EI50), gam = as.numeric(best_model$gam))
      )

      tryCatch({
        fit <- nls(model_formula, start = start_list, data = dataset)

        # Calculate predicted y values manually
        x_vals <- seq(min(dataset[[index]]), max(dataset[[index]]), length.out = 100)
        y_vals <- best_model$E0 - best_model$Emax * x_vals ^ best_model$gam /
          (best_model$EI50 ^ best_model$gam + x_vals ^ best_model$gam)

        # Custom x-axis label: Check if x_label is a list or a single value
        if (is.null(x_label)) {
          c_x_label <- index  # Default: set c_x_label as the current x_column name
        } else if (is.list(x_label) && !is.null(x_label[[index]])) {
          c_x_label <- x_label[[index]]  # Use custom label for the current x_column
        } else {
          c_x_label <- x_label  # Use global custom label if provided
        }

        # Check if log_scale_x is TRUE for the current index
        log_scale_x_for_index <- ifelse(is.list(log_scale_x) && index %in% names(log_scale_x), log_scale_x[[index]], FALSE)

        if (log_scale_x_for_index) {
          # Dynamically adjust zero values based on the lowest non-zero value
          non_zero_values <- dataset[[index]][dataset[[index]] > 0]
          if (length(non_zero_values) > 0) {
            min_non_zero <- min(non_zero_values)
            adjusted_zero_value <- 10^(floor(log10(min_non_zero)) - 1)

            zero_count <- sum(dataset[[index]] == 0)
            if (zero_count > 0) {
              message(paste("Adjusted", zero_count, "zero values in", index,
                            "to", adjusted_zero_value,
                            "(1 log10 lower than the smallest non-zero value)."))
              dataset[[index]][dataset[[index]] == 0] <- adjusted_zero_value
            }
          } else {
            stop("No non-zero values found in the dataset for index: ", index)
          }

          # Apply log scale for x-axis
          c_x_label <- paste("Log10(", index, ")", sep = "")

          # Adjust the ylim by adding 1 to the maximum value if it's not provided
          if (is.null(ylim)) {
            ylim <- c(min(dataset[[y_column]]), max(dataset[[y_column]]) + 1)  # Extend the y-axis by 0.5
          } else {
            ylim[2] <- ylim[2] + 1  # Extend the upper limit of provided ylim by 1
          }

          # Plot data with log scale for x-axis
          plot(dataset[[index]], dataset[[y_column]],
               main = c_plot_title,
               xlab = c_x_label, ylab = c_y_label,
               pch = 19, col = dot_color,
               xlim = xlim, ylim = ylim,
               log = "x", xaxt = "n", # Disable default x-axis
               cex.main = title_cex,
               cex.lab  = label_cex,
               cex.axis = axis_cex)

          # Custom ticks for the log scale (powers of 10)
          axis(1, at = 10^seq(floor(log10(min(dataset[[index]]))),
                              ceiling(log10(max(dataset[[index]])))), labels = TRUE, cex.axis = axis_cex)
        } else {
          # If no log scale is specified, plot as normal
          if (is.null(ylim)) {
            ylim <- c(min(dataset[[y_column]]), max(dataset[[y_column]]) + 1)  # Extend the y-axis by 0.5
          } else {
            ylim[2] <- ylim[2] + 1  # Extend the upper limit of provided ylim by 1
          }

          plot(dataset[[index]], dataset[[y_column]],
               main = c_plot_title,
               xlab = c_x_label, ylab = c_y_label,
               pch = 19, col = dot_color,
               xlim = xlim, ylim = ylim,
               cex.main = title_cex,
               cex.lab  = label_cex,
               cex.axis = axis_cex)
        }

        # Overlay the fitted model curve
        lines(x_vals, y_vals, col = fit_color, lwd = 2)
        # Extract and convert R_squared percentage back to numeric if it's in character format
        R2_numeric <- ifelse(!is.na(best_model$R_squared) && grepl("%", best_model$R_squared),
                             as.numeric(gsub("%", "", best_model$R_squared)) / 100,
                             NA)

        R2_val <- ifelse(!is.na(R2_numeric), paste0(round(R2_numeric * 100, 2), "%"), "NA")

        y_max <- ylim[2]  # Get the upper limit of y-axis
        y_position <- y_max * 0.75  # Position the text at 85% of the y-axis max value (you can adjust this value)
        x_position <- max(dataset[[index]]) * 0.1  # Position the text at 10% of the max x-axis value (you can adjust this too)

        # Model details on the plot
        model_details <- paste(
          "AIC: ", round(best_model$AIC, 2), "\n",
          "R\u00B2: ", R2_val, "\n",
          "E0: ", round(best_model$E0, 4), "\n",
          "Emax: ", round(best_model$Emax, 4), "\n",
          "Gam: ", round(best_model$gam, 4), "\n",
          "EI50: ", round(best_model$EI50, 4)
        )

        # Display model details on the plot
        text(x = x_position, y = y_position, labels = model_details, pos = 4, cex = detail_cex, col = "black")


        # Store the plot in the list
        plots[[index]] <- recordPlot()

      }, error = function(e) {
        message("Plot skipped for index ", index, ": ", e)
      })
    }
  }

  return(list(All_Model_Results = all_results, Best_Models = best_models,Plots = plots))
}

