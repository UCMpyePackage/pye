#' @title Center and Scale Continuous Regressors for pye Modeling
#'
#' @description Scales and centers specified numeric variables in a data
#'   frame (mean 0, SD 1). It automatically excludes the target variable
#'   (\code{y}) and variables identified as binary (dummy), ensuring only
#'   continuous regressors are standardized. This prepares data for models
#'   where standardization is necessary, such as the pye optimization.
#'
#' @param df A \code{data.frame} or \code{tibble} containing all variables.
#' @param X A \code{character} vector of column names considered as potential
#'   regressors. Only numeric variables from this list will be scaled.
#' @param y A \code{character} string: the name of the target variable.
#'   This variable will not be scaled.
#'
#' @return A \code{list} containing two elements:
#'   \item{df_scaled}{A \code{data.frame} with the selected numeric
#'     regressors scaled (centered to mean 0 and scaled to SD 1).}
#'   \item{scaling_params}{A \code{list} containing:
#'     \itemize{
#'       \item \code{original_mu}: Means of the variables that were scaled.
#'       \item \code{original_stdev}: Standard deviations of the variables.
#'     }
#'     These parameters enable consistent scaling of new data.}
#'
#' @details
#' The function first identifies potential regressors from the \code{X} argument.
#' It then excludes variables that are the target (\code{y}), are not
#' numeric, or are identified as binary (having exactly two unique values).
#' Only the remaining continuous numeric variables are scaled.
#' The scaling parameters are returned to allow for the consistent
#' transformation of future datasets (e.g., test or new data).
#'
#' @examples
#' library(pye)
#'
#' # Simulate the dataframe
#' sim_data <- create_sample_with_covariates(
#'   rows_train = 100, rows_test = 50, cols = 40, cols_cov = 12, max_rho = 0.3, seed = 1)
#' df <- sim_data$df
#' X_cols_name <- sim_data$X
#' y_col_name <- sim_data$y
#' C_cols_name <- sim_data$C
#'
#' # Now call the scaling function
#' scaled_output <- scaling_df_for_pye(
#'   df = df,
#'   X = c(X_cols_name, C_cols_name),
#'   y = y_col_name
#' )
#'
#' print(head(scaled_output$df_scaled))
#' print(scaled_output$scaling_params$original_mu)
#' print(scaled_output$scaling_params$original_stdev)
#'
#' @importFrom stats sd na.omit
#' @importFrom methods is
#' @export
scaling_df_for_pye <- function(df, X, y) {

  # --- Input Validation ---
  if (!is.data.frame(df) || nrow(df) == 0) stop("Input data frame 'df' must be a non-empty data frame.")
  if (!is.character(X) || length(X) == 0) {stop("Input 'X' must be a non-empty character vector of column names.")}
  if (!is.character(y) || length(y) != 1) {stop("Input 'y' must be a single character string representing the target column name.")}

  # Check if all specified X and y columns exist in df
  missing_X <- setdiff(X, names(df))
  if (length(missing_X) > 0) {
    stop("The following 'X' columns are not found in 'df': ",
         paste(missing_X, collapse = ", "))
  }
  if (!(y %in% names(df))) {
    stop(paste0("The target variable '", y, "' is not found in 'df'."))
  }

  # --- Identify columns not to scale ---
  # 1. The target variable 'y'
  # 2. Columns in 'df' that are NOT in the provided 'X' regressors list
  # 3. Columns within 'X' that are identified as binary (have exactly 2 unique values)
  #    Ensure this check only happens for columns actually present in 'df' and 'X'

  # Columns that are truly potential regressors (in X and in df)
  potential_regressors_in_df <- intersect(X, names(df))

  # Identify binary variables among potential regressors
  binary_vars <- character(0)
  if (length(potential_regressors_in_df) > 0) {
    # Use sapply with tryCatch for robust handling of potential errors with unique()
    # e.g., if a column is factor with NA, unique() might return 3 values.
    # This check specifically targets numeric or simple character/factor with 2 levels.
    is_binary <- sapply(df[, potential_regressors_in_df, drop = FALSE], function(col) {
      length(unique(stats::na.omit(col))) == 2
    })
    binary_vars <- names(is_binary[is_binary])
  }

  # All columns that should explicitly NOT be scaled
  to_not_scale <- unique(c(y, setdiff(names(df), X), binary_vars))

  # Identify numeric columns among the remaining candidates for scaling
  # These are the columns in df, that are in X, are not y, and are not binary, AND are numeric.

  # First, find candidate columns that are in X and not in 'to_not_scale'
  candidate_for_scaling_names <- setdiff(names(df), to_not_scale)

  # Then, filter these candidates to only include numeric columns
  vars_to_scale <- names(df[, candidate_for_scaling_names, drop = FALSE])[sapply(df[, candidate_for_scaling_names, drop = FALSE], is.numeric)]

  if (length(vars_to_scale) == 0) {
    warning("No numeric variables were identified for scaling based on the provided inputs. Returning original df.")
    return(list(df_scaled = df, original_mu = numeric(0), original_stdev = numeric(0)))
  }

  # --- Perform Scaling ---
  original_var_mean <- colMeans(df[, vars_to_scale, drop = FALSE], na.rm = TRUE)
  original_var_stdev <- sapply(df[, vars_to_scale, drop = FALSE], stats::sd, na.rm = TRUE)

  df_scaled <- df # Create a copy to modify
  # Apply scale function to selected numeric columns
  df_scaled[, vars_to_scale] <- as.data.frame(
    scale(df[, vars_to_scale, drop = FALSE], center = TRUE, scale = TRUE)
  )

  return(list(df_scaled = df_scaled,
              scaling_params = list(original_mu = original_var_mean,
                                    original_stdev = original_var_stdev)))
}