#-- Gradient of the function: phi(vtheta):=-(logL(vtheta) / T) --#
fn_gr <- function(v_theta, l_data) {
     # retrieve data
     m_y  <- l_data[["m_y" ]]
     m_ys <- l_data[["m_ys"]]
     l_x  <- l_data[["l_x" ]]
     m_W  <- l_data[["m_W" ]]

     # sample size
     N <- nrow(m_y)
     TT <- ncol(m_y)
     K <- length(l_x)

     v_psi   <- v_theta[1:N] # (N,)
     v_beta  <- v_theta[(N + 1):(N + (K * N))]; # (KN,) x 1
     v_sgmsq <- v_theta[(N + (K * N) + 1):(N + (N * K) + N)] # (N,)

     # reshape beta
     m_beta_tr <- matrix(v_beta, K, N) # (K,N)
     m_beta <- t(m_beta_tr) # (N,K)

     v_sgm4h <- v_sgmsq^2
     v_sgm6h <- v_sgmsq^3

     # compute residuals
     m_beta_times_x <- matrix(0, N, TT) # 0s needed for recursive sum; (N,T) with generic element \vbeta_{i}'\vx_{it}
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          v_beta_k <- m_beta[, k] # (N,)
          m_beta_times_x_k <- m_x_k * v_beta_k #!! (N,T) "times" (N,1) = (N,T)
          m_beta_times_x <- m_beta_times_x + m_beta_times_x_k
     }
     m_psi_times_ys <- m_ys * v_psi #!! (N,T) "times" (N,1) = (N,T)
     m_eps <- m_y - m_psi_times_ys - m_beta_times_x # N x T
     v_ssr <- rowSums(m_eps^2) # (N,) sum of squared residuals: sum_t (eps_it)^2 !!crossprod()??

     # standardize
     sssr <- sum(v_ssr / v_sgmsq) # sum of squared standardized-residuals: sum_i sum_t (eps_it / sgm_i)^2

     m_Psi <- diag(v_psi, nrow = N) # (N,N)
     m_A <- diag(N) - (m_Psi %*% m_W)

     det_mA <- det(m_A)
     if (det_mA <= 0) {
          stop("`I - Psi W` is not invertible. Check that the parameter space is defined correctly.")
     }

     # first derivative
     ## v_dphi_dvpsi
     m_Q <- m_W %*% solve(m_A) #!! IS THERE A BETTER WAY TO WRITE THIS LINE?
     v_dphi_dvpsi <- diag(m_Q) - (rowSums(m_ys * m_eps) / v_sgmsq / TT) #!! stopifnot(dim(m_Q) > 1)

     ## v_dphi_dvbeta
     m_X_times_eps_divided_sgmsq <- matrix(NA_real_, N, K)
     m_eps_divided_sgmsq <- m_eps / v_sgmsq # (N,T)
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          m_x_k_times_eps_divided_sgmsq <- m_x_k * m_eps_divided_sgmsq # (N,T)
          m_X_times_eps_divided_sgmsq[, k] <- rowSums(m_x_k_times_eps_divided_sgmsq) # (N,)
     }
     m_X_times_eps_divided_sgmsq_tr <- t(m_X_times_eps_divided_sgmsq) # (K,N)
     v_dphi_dvbeta <- -c(m_X_times_eps_divided_sgmsq_tr) / TT # (NK,)

     ## v_dphi_dvsgmsq
     v_dphi_dvsgmsq <- (1 / 2 / v_sgmsq) - (v_ssr / v_sgm4h / TT / 2)
     v_gradient <- c(v_dphi_dvpsi,
                     v_dphi_dvbeta,
                     v_dphi_dvsgmsq)

     # return
     v_gradient
}
