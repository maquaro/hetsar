#' Heterogeneous SAR estimation
#'
#' Estimation of heterogeneous spatial autoregressive models
#'
#' @param data A `data.frame`, in tidy format. The panel _must_ be balanced.
#' @param W Spatial weight `matrix`. Rows and columns should be ordered such that colnames(W)==rowanems(W)==unique(df$id).
#' @param indices Character vector of length 2 with the names of the columns in `data` containing the unit-identifier (`indices[[1]]`) and the time-identifier (`indices[[2]]`).
#' @param dependent Name of the dependent variable.
#' @param explanatory Character vector with the names of the explanatory variables, excluding: (a) the spatial lag of the dependent variable `Wy`, included by default; (b) the intercept; (c) space/time-lags. See Details.
#' @param nocons Logical: Shall the constant term (i.e. the intercept) be excluded?
#' @param time_lags Character vector with the names of the time-lagged variables.
#' @param space_lags Character vector with the names of the space-lagged variables.
#' @param time_space_lags Character vector with the names of the time- and space-lagged variables.
#' @param verbose Character vector. Verbosely print some extra information relative to the model specification (`verbose="model"`), to the bounds (`verbose="bounds"`), to the optimisation procedure (`verbose="optimisation"`), or to all the above (`verbose="all"). Default is silent (`verbose="silent"`).
#' @param bounds (Advanced) Named-`list` specifying the parameter space, see "Bounds".
#'
#' @return
#' TODO An object (`list`) of class `hetsar` with the following elements:
#' \item{coefficients}{A matrix of estimated coefficients. !! COLNAMES}
#' \item{Std.err.sandwich}{A matrix with standard errors computed using the sandwich formula.}
#' \item{Std.err.standard}{A matrix with standard errors computed using the standard formula.}
#' \item{N}{The number of spatial units.}
#' \item{T}{The number of time periods.}
#' \item{NT}{The number of observations.}
#' 
#' @section Model specification:
#' TODO
#' 
#' @section Bounds:
#' TODO
#'
#' @seealso
#' TODO See also \code{\link[hetsar]{summary.hetsar}} to see the results.
#'
#' @author
#' Michele Aquaro
#'
#' @references
#' M. Aquaro, N. Bailey and M. H. Pesaran (2021). "Estimation and inference for spatial models with heterogeneous coefficients: An application to US house prices". Journal of Applied Econometrics 36(1): 18-44. doi:https://doi.org/10.1002/jae.2792
#'
#' @examples
#'
#' data(hetsarDataDemo)
#' df_data <- hetsarDataDemo[["data"]]
#' m_C <- hetsarDataDemo[["network_matrix"]]
#' m_W <- m_C / rowSums(m_C) # row normalise
#'
#' out <- hetsar(
#'   data = df_data, 
#'   W = m_W,
#'   indices = c("id", "time"),
#'   dependent = "y",
#'   explanatory = "x",
#'   space_lags = "x")
#'
#' summary(out)
#' summary(out, MG = TRUE)
#'
#' @export
hetsar <- function(
     data, 
     W, 
     indices, # can [[1]] be a character? Can [[2]] be a date?
     dependent,
     explanatory = NULL,
     nocons = FALSE,
     time_lags = NULL,
     space_lags = NULL,
     time_space_lags = NULL,
     verbose = "silent",
     bounds = NULL) {

     # -------------------------------------------------------------------------- #
     # rename
     df_data <- data
     m_W <- W
     #!! rm(list(data, m_W)) ?

     # ========================================================================== #
     # CHECKS

     # -------------------------------------------------------------------------- #
     # checks type/class

     if (!any(class(df_data) == "data.frame")) {
          stop("`data` must be an object of class 'data.frame'")
     }

     if (!any(class(m_W) == "matrix")) {
          stop("`W` must be an object of class 'matrix'")
     }

     if (!typeof(dependent) == "character") {
          stop("`dependent` must be an object of type 'character'")
     }

     if (!is.null(explanatory)) {
          if (!typeof(explanatory) == "character") {
               stop("`explanatory` must be an object of type 'character'")
          }
     }

     if (!typeof(nocons) == "logical") {
          stop("`nocons` must be an object of type 'logical'")
     }

     if (!is.null(time_lags)) {
          if (!typeof(time_lags) == "character") {
               stop("`time_lags` must be an object of type 'character'")
          }
     }

     if (!is.null(space_lags)) {
          if (!typeof(space_lags) == "character") {
               stop("`space_lags` must be an object of type 'character'")
          }
     }

     if (!is.null(time_space_lags)) {
          if (!typeof(time_space_lags) == "character") {
               stop("`time_space_lags` must be an object of type 'character'")
          }
     }

     if (!typeof(verbose) == "character") {
          stop("`verbose` must be an object of type 'character'")
     }

     if (!is.null(bounds)) {
          if (!typeof(bounds) == "list") {
               stop("`bounds` must be an object of type 'list'")
          }
     }

     # -------------------------------------------------------------------------- #
     # other checks

     if (length(dependent) > 1) {
          stop("`dependent` _must_ be a character vector of length 1 (i.e. only one outcome variable is allowed)")
     }

     if (length(indices) != 2) {
          stop("One of the two indices is missing (`indices` _must_ be a character vector of length 2)")
     }

     # -------------------------------------------------------------------------- #
     # checks on `df_data`

     # 1. are all variables to be used in data?
     v_variables <- c(indices, dependent, explanatory, time_lags, space_lags, time_space_lags)
     v_bool <- v_variables %in% names(df_data)
     if (!all(v_bool)) {
          # location of the unmatched column name(s)
          v_ind <- which(!v_bool)
          stop(sprintf("%s not found in `data`", paste(v_variables[v_ind], collapse = ", ")))
     }
     df_data <- df_data[, v_variables] # copy

     # 2. remove NAs
     NT_old <- nrow(df_data)
     df_data <- stats::na.omit(df_data)
     NT_new <- nrow(df_data)
     if (NT_new < NT_old) {
          cat(sprintf("  Removed %d row(s) with NAs\n", NT_old - NT_new))
     }

     # 3. balanced?
     v_NT_idx_i <- df_data[[indices[[1L]]]] # (NT,)
     v_NT_idx_t <- df_data[[indices[[2L]]]]
     N. <- length(unique(v_NT_idx_i))
     T. <- length(unique(v_NT_idx_t))
     NT <- N. * T.
     if (nrow(df_data) != NT) {
          stop("`data` is not balanced (i.e. nrow(df_data)!=NT)")
     }

     # 4. sort
     df_data <- df_data[order(v_NT_idx_t, v_NT_idx_i), ]

     # 5. no "L." or "W." or "L.W." in the names
     v_colnames <- names(df_data)
     v_ind <- grep("L.", substr(v_colnames, 1, 2), fixed = TRUE)
     if (length(v_ind) > 0) {
          stop(sprintf("Rename column(s) %s in `data` (pls use something different than `L.foo`, see Documentation)", paste(v_colnames[v_ind], collapse = ", ")))
     }
     v_ind <- grep("W.", substr(v_colnames, 1, 2), fixed = TRUE)
     if (length(v_ind) > 0) {
          stop(sprintf("Rename column(s) %s in `data` (pls use something different than `W.foo`, see Documentation)", paste(v_colnames[v_ind], collapse = ", ")))
     }

     # -------------------------------------------------------------------------- #
     # checks on `m_W`
     
     if (nrow(m_W) != ncol(m_W)) {
          stop("`W` must be a square matrix")
     }

     if (any(c(m_W) < 0)) {
          stop("negative weights in `m_W` are not allowed")
     }

     if (!all(rowSums(m_W) > 0)) {
          stop("each spatial unit in `m_W` must have at least 1 neighbour")
     }

     if (nrow(m_W) != N.) {
          stop("`df_data` and `m_W` must have the same number of spatial units N")
     }

     # -------------------------------------------------------------------------- #
     # checks `l_x` _must_ have at least one element
     v_bool <- c(
          nocons,
          is.null(explanatory),
          is.null(time_lags),
          is.null(space_lags),
          is.null(time_space_lags))
     if (all(v_bool)) {
          stop("'y = Wy + e' not coded. Please add at least one variable on the rhs of the equation")
     }

     # -------------------------------------------------------------------------- #
     # checks on bounds

     # construct the full list of names that can appear in bounds
     v_name_bounds_universe <- c(paste0("W.", dependent), "sigmasq")
     if (!nocons) {
          v_name_bounds_universe <- c(v_name_bounds_universe, "cons")
     }
     if (!is.null(explanatory)) {
          v_name_bounds_universe <- c(v_name_bounds_universe, explanatory)
     }
     if (!is.null(time_lags)) {
          v_name_bounds_universe <- c(v_name_bounds_universe, paste0("L.", time_lags))
     }
     if (!is.null(space_lags)) {
          v_name_bounds_universe <- c(v_name_bounds_universe, paste0("W.", space_lags))
     }
     if (!is.null(time_space_lags)) {
          v_name_bounds_universe <- c(v_name_bounds_universe, paste0("L.W.", time_space_lags))
     }
     # (do I need to worry about duplicates? Perhaps not)
     v_name_bounds_universe <- unique(v_name_bounds_universe)
     if (!is.null(bounds)) {
          v_bool <- names(bounds) %in% v_name_bounds_universe
          if (!all(v_bool)) {
               stop(sprintf("%s %s",
                    "Some names in `bounds` seem to be mispelled (i.e. `names(bounds)`).\n  Allowed names in `bounds` are: ",
                    paste(v_name_bounds_universe, collapse = ", ")))
          }
     }

     # ========================================================================== #
     # LONG TO WIDE
     
     # -------------------------------------------------------------------------- #
     # explanatory variables
     l_x <- list()

     ## constant
     if (!nocons) {
          ## remove t=1 if there are time-lags
          if (length(time_lags) > 0 || length(time_space_lags) > 0) {
               l_x[["cons"]] <- matrix(1, N., T. - 1)
          } else {
               l_x[["cons"]] <- matrix(1, N., T.)
          }
     }

     ## non-lagged explanatory variables
     if (!is.null(explanatory)) {
          for (x in explanatory) {
               m_x <- matrix(df_data[[x]], N., T.)

               ## remove t=1 if there are time-lags
               if (length(time_lags) > 0 || length(time_space_lags) > 0) {
                    l_x[[x]] <- m_x[, 2:T.]
               } else {
                    l_x[[x]] <- m_x
               }
          }
     }

     ## time-lagged explanatory variables
     if (!is.null(time_lags)) {
          for (x in time_lags) {
               m_x <- matrix(df_data[[x]], N., T.)
               l_x[[paste0("L.", x)]] <- m_x[, 1:(T. - 1)]
          }
     }

     ## space-lagged explanatory variables
     if (!is.null(space_lags)) {
          for (x in space_lags) {
               m_x <- matrix(df_data[[x]], N., T.)
               m_Wx <- m_W %*% m_x

               ## remove t=1 if there are time-lags
               if (length(time_lags) > 0 || length(time_space_lags) > 0) {
                    l_x[[paste0("W.", x)]] <- m_Wx[, 2:T.]
               } else {
                    l_x[[paste0("W.", x)]] <- m_Wx
               }
          }
     }

     ## time-space-lagged explanatory variables
     if (!is.null(time_space_lags)) {
          for (x in time_space_lags) {
               m_x <- matrix(df_data[[x]], N., T.)
               m_Wx <- m_W %*% m_x
               l_x[[paste0("L.W.", x)]] <- m_Wx[, 1:(T. - 1)]
          }
     }

     # -------------------------------------------------------------------------- #
     # dependent variable
     m_y <- matrix(df_data[[dependent]], N., T.) 

     # remove t=1 if there are time-lags
     if (length(time_lags) > 0 || length(time_space_lags) > 0) {
          m_y <- m_y[, 2:T.]
     }

     # -------------------------------------------------------------------------- #
     # pack all data into a list
     l_data <- list()
     l_data[["m_y"]] <- m_y
     l_data[["l_x"]] <- l_x
     l_data[["m_W"]] <- m_W
     #browser()

     # ========================================================================== #
     # bounds

     K <- length(l_x)
     N <- N.

     # character vector with parameter names
     v_name_parameter <- c(
          paste0("W.", dependent),
          names(l_x),
          "sigmasq")

     # if `bounds` is NULL assign it to a named-list
     if (is.null(bounds)) {
          bounds <- vector("list", length(v_name_parameter)) # (K+2,)
          names(bounds) <- v_name_parameter
     }

     bnd <- 0.995 #cc#

     # default lower/upper-bounds for W.y
     if (is.null(bounds[[paste0("W.", dependent)]])) {
          bounds[[paste0("W.", dependent)]] <- c(-bnd, bnd)
     }

     # default lower/upper-bounds for sigmasq
     if (is.null(bounds[["sigmasq"]])) {
          bounds[["sigmasq"]] <- c(0.01, Inf)
     }

     # default lower/upper-bounds for L.y
     if ( (dependent %in% time_lags) &&  is.null(bounds[[paste0("L.", dependent)]]) ) {
          bounds[[paste0("L.", dependent)]] <- c(-bnd, bnd)
     }

     # default lower/upper-bounds for L.W.y
     if ( (dependent %in% time_lags) &&  is.null(bounds[[paste0("L.W.", dependent)]]) ) {
          bounds[[paste0("L.W.", dependent)]] <- c(-bnd, bnd)
     }

     # do the rest of the variables (e.g. `x`, `W.x`)
     for (np in v_name_parameter) { # mnemonic: Name Parameter
          if (is.null(bounds[[np]])) {
               # default lower/upper-bounds for x
               bounds[[np]] <- c(-Inf, Inf)
          } else {
               next
          }
     }
     #browser()

     ## from "list" to lower/upper-bound "matrix"
     m_lb_theta <- matrix(NA_real_, N, K + 2)
     m_ub_theta <- matrix(NA_real_, N, K + 2)
     colnames(m_lb_theta) <- v_name_parameter
     colnames(m_ub_theta) <- v_name_parameter
     for (np in v_name_parameter) { # mnemonic: Name Parameter
          vec <- bounds[[np]] # (2,)
          m_lb_theta[, np] <- rep(vec[[1]], N)
          m_ub_theta[, np] <- rep(vec[[2]], N)
     }
     #print(m_lb_theta)
     #print(m_ub_theta)

     ## reshape m_lb_beta
     m_lb_beta <- m_lb_theta[, names(l_x)] # (N,K)
     m_ub_beta <- m_ub_theta[, names(l_x)] # (N,K)
     m_lb_beta_tr <- t(m_lb_beta) # (K,N)
     m_ub_beta_tr <- t(m_ub_beta) # (K,N)
     v_lb_beta <- c(m_lb_beta_tr) # (NK,)
     v_ub_beta <- c(m_ub_beta_tr) # (NK,)

     ## fill v_lb and v_ub
     v_lb <- vector("numeric", N * (K + 2))
     v_ub <- vector("numeric", N * (K + 2))

     v_lb[1:N] <- m_lb_theta[, paste0("W.", dependent)] # (N,)
     v_ub[1:N] <- m_ub_theta[, paste0("W.", dependent)] # (N,)

     v_lb[(N + 1):(N + (K * N))] <- v_lb_beta  # (KN,)
     v_ub[(N + 1):(N + (K * N))] <- v_ub_beta  # (KN,)

     v_lb[(N + (K * N) + 1):(N + (N * K) + N)] <- m_lb_theta[, "sigmasq"] # (N,)
     v_ub[(N + (K * N) + 1):(N + (N * K) + N)] <- m_ub_theta[, "sigmasq"] # (N,)

     l_bounds <- list(v_lb, v_ub)
     #browser()

     # ========================================================================== #
     # VERBOSE

     cat("SAR with heterogeneous coefficients\n\n")

     # formula
     if (verbose == "all" || verbose == "model") {

          # formula
          Kp2 <- length(v_name_parameter)
          cat("Formula (short-hand notation):\n")
          cat(sprintf("  %s = %s + \u03b5\n\n", dependent, paste(v_name_parameter[1:(Kp2 - 1)], collapse = " + ")))
     }

     # bounds
     if (verbose == "all" || verbose == "bounds") {
          cat("Bounds (e.g. parameter space):\n")
          df_temp <- data.frame(variable = v_name_parameter)
          for (np in v_name_parameter) {
               vec <- bounds[[np]] # (2,)
               df_temp[df_temp[["variable"]] == np, "lower_bounds"] <- vec[[1]]
               df_temp[df_temp[["variable"]] == np, "upper_bounds"] <- vec[[2]]
          }
          print(df_temp)
     }

     # ========================================================================== #
     # HSAR
     l_out <- fn_hsar(l_data = l_data, l_bounds = l_bounds)
     m_theta <- l_out[["theta"]]
     m_se_sandwich <- sqrt(l_out[["sandwich"]])
     m_se_standard <- sqrt(l_out[["standard"]])

     # residuals
     m_resid_rform <- l_out[["resid_rform"]]
     m_resid_naive <- l_out[["resid_naive"]]
     if (length(time_lags) > 0 || length(time_space_lags) > 0) {
          # remove t=1 if there are time-lags
          df_resid <- df_data[(N + 1):nrow(df_data), indices] # (NT,2)
     } else {
          df_resid <- df_data[, indices] # (NT,2)
     }
     df_resid[["reduced_form"]] <- c(m_resid_rform)
     df_resid[["naive"]]        <- c(m_resid_naive)

     # report verbosely information on optimisation
     if (verbose == "all" || verbose == "optimisation") {
          cat("\n  Optimisation (see `?optim`):\n")
          cat(sprintf("%s: %d\n", "  - exit flag", l_out[["exitflag"]]))
          cat(sprintf("%s: %f\n", "  - value of the objective function (-1/T)", l_out[["fval"]]))
          cat(sprintf("%s: %s\n", "  - any other message from the solver", l_out[["optim.msg"]]))
     }

     # ========================================================================== #
     # prepare output

     #!! TODO shall I leave these as matrices or shall I transform them into data-frames? 
     #!! What would be the advantage of using data-frames?
     #!! - if I have to re-shuffle columns

     v_id <- unique(df_data[[indices[[1L]]]]) # (N,)
     rownames(m_theta)       <- v_id
     rownames(m_se_sandwich) <- v_id
     rownames(m_se_standard) <- v_id

     colnames(m_theta)       <- v_name_parameter
     colnames(m_se_sandwich) <- v_name_parameter
     colnames(m_se_standard) <- v_name_parameter

     l_return <- list()
     l_return[["coefficients"]] <- m_theta
     l_return[["Std.err.sandwich"]] <- m_se_sandwich
     l_return[["Std.err.standard"]] <- m_se_standard
     l_return[["resid"]] <- df_resid
     l_return[["T"]] <- ncol(m_y)
     l_return[["N"]] <- nrow(m_y)
     l_return[["NT"]] <- nrow(m_y) * ncol(m_y)
     l_return[["fval"]] <- l_out[["fval"]]
     l_return[["exitflag"]] <- l_out[["exitflag"]]
     l_return[["optim.msg"]] <- l_out[["optim.msg"]]

     # add class to the S3 object
     class(l_return) <- "hetsar"

     l_return
}
