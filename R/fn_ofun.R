#-- Objective function: phi(vtheta):=-(logL(vtheta) / T) --#
fn_ofun <- function(v_theta, l_data) {
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
     m_beta_times_x <- matrix(0, N, TT) # 0 needed for recursive sum; (N x T) with generic element \vbeta_{i}'\vx_{it}
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          v_beta_k <- m_beta[, k] # (N,)
          m_beta_times_x_k <- m_x_k * v_beta_k #!!
          #browser()
          #print(dim(m_beta_times_x))
          #print(dim(m_beta_times_x_k))
          m_beta_times_x <- m_beta_times_x + m_beta_times_x_k
     }
     m_psi_times_ys <- m_ys * v_psi #!! N x T
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
     constant <- log(2 * pi) * N / 2
     first_part <- -log(det_mA)
     secon_part <- sum(log(v_sgmsq)) / 2
     third_part <- sssr / TT / 2
     ofun <- constant + first_part + secon_part + third_part
     #print(ofun)
     #print(first_part)
     #print(secon_part)
     #print(third_part)
     #print(v_sgmsq)
     #print(sssr)
     #print(v_ssr)
     #print(v_sgmsq)

     # return
     ofun
}
