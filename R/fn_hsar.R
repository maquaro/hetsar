#-- wrapper: write the estimator as a function of the dataset, only --#
fn_hsar <- function(l_data, l_bounds) {
     # -------------------------------------------------------------------------- #
     # l_data is a list containing:
     # m_y: (N,T) matrix
     # l_x: list containing K matrices of covariates of order (N,T) (it may or may note include the intercept)
     # m_W: (N,N) matrix

     # -------------------------------------------------------------------------- #
     N  <- nrow(l_data[["m_y"]])
     TT <- ncol(l_data[["m_y"]])
     K <- length(l_data[["l_x"]])

     # add y*:=Wy
     l_data[["m_ys"]] <- l_data[["m_W"]] %*% l_data[["m_y"]] # (N,T)
     

     # -------------------------------------------------------------------------- #
     # constructs v_theta_ini
     l_theta_ini <- NULL #cc#
     v_theta_ini <- vector("numeric", N * (K + 2))
     if (is.null(l_theta_ini)) {
          v_theta_ini[1:(N + (K * N))] <- 0
          v_theta_ini[(N + (K * N) + 1):(N + (N * K) + N)] <- 1
     } else {
          v_psi_ini   <- l_theta_ini[["psi"  ]]
          m_beta_ini  <- l_theta_ini[["beta" ]]
          v_sgmsq_ini <- l_theta_ini[["sqmsq"]]

          ## reshape m_beta_ini
          m_beta_ini_tr <- t(m_beta_ini) # (K,N)
          v_beta_ini <- c(m_beta_ini_tr) # (NK,)

          ## fill v_theta_ini
          v_theta_ini[1:N]                                 <- v_psi_ini   # (N,)
          v_theta_ini[(N + 1):(N + (K * N))]               <- v_beta_ini  # (KN,)
          v_theta_ini[(N + (K * N) + 1):(N + (N * K) + N)] <- v_sgmsq_ini # (N,)
     }

     # -------------------------------------------------------------------------- #
     l_out <- stats::optim(
          par     = v_theta_ini, 
          fn      = fn_ofun,
          gr      = fn_gr,
          l_data  = l_data,
          method  = "L-BFGS-B",
          lower   = l_bounds[[1]], # (N(K+2),)
          upper   = l_bounds[[2]], # (N(K+2),)
          control = list(maxit = 1e6))
     v_theta  <- l_out[["par"]]
     #print(l_out)

     v_psi   <- v_theta[1:N] # (N,)
     v_beta  <- v_theta[(N + 1):(N + (K * N))]; # (KN,)
     v_sgmsq <- v_theta[(N + (K * N) + 1):(N + (N * K) + N)] # (N,)
     #print(range(v_sgmsq))

     # reshape beta
     m_beta_tr <- matrix(v_beta, K, N) # (K,N)
     m_beta <- t(m_beta_tr) # (N,K)

     # pack all estimates into a matrix with names
     m_theta <- matrix(NA_real_, N, K + 2)
     m_theta[, 1]         <- v_psi
     m_theta[, 2:(K + 1)] <- m_beta
     m_theta[, (K + 2)]   <- v_sgmsq
     #browser()
     #print(m_theta)

     # -------------------------------------------------------------------------- #
     # residuals

     m_y  <- l_data[["m_y" ]]
     m_ys <- l_data[["m_ys"]]
     l_x  <- l_data[["l_x" ]]
     m_W  <- l_data[["m_W" ]]
     N <- nrow(m_y)
     TT <- ncol(m_y)
     K <- length(l_x)

     m_beta_times_x <- matrix(0, N, TT) # 0 needed for recursive sum; (N x T) with generic element \vbeta_{i}'\vx_{it}
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          v_beta_k <- m_beta[, k] # (N,)
          m_beta_times_x_k <- m_x_k * v_beta_k #!!
          m_beta_times_x <- m_beta_times_x + m_beta_times_x_k
     }
     m_psi_times_ys <- m_ys * v_psi #!! N x T
     m_yhat_naive <- m_psi_times_ys + m_beta_times_x
     m_resid_naive <- m_y - m_yhat_naive # (N,T)

     m_Psi <- diag(v_psi, nrow = N) # (N,N)
     m_A <- diag(N) - (m_Psi %*% m_W)
     m_yhat_rform <- solve(m_A, m_beta_times_x)
     m_resid_rform <- m_y - m_yhat_rform # (N,T)

     # -------------------------------------------------------------------------- #
     # variance
     l_variance <- fn_var_ml(m_theta, l_data) # (2,)

     # -------------------------------------------------------------------------- #
     # return
     l_return <- list()
     l_return[["theta"]] <- m_theta # (3,)

     l_return[["standard"]] <- l_variance[["standard"]]
     l_return[["sandwich"]] <- l_variance[["sandwich"]]

     l_return[["resid_rform"]] <- m_resid_rform
     l_return[["resid_naive"]] <- m_resid_naive

     l_return[["fval"]] <- l_out[["value"]]
     l_return[["exitflag"]] <- l_out[["convergence"]]
     l_return[["optim.msg"]] <- l_out[["message"]]

     l_return
}
