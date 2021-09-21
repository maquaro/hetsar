#-- \hat{\var}(\mTheta), where m_theta = (v_psi|m_beta|v_sgmsq) --#
fn_var_ml <- function(m_theta, l_data) {
     # retrieve data
     m_y  <- l_data[["m_y" ]]
     m_ys <- l_data[["m_ys"]]
     l_x  <- l_data[["l_x" ]]
     m_W  <- l_data[["m_W" ]]

     # sample size
     N  <- nrow(m_y)
     TT <- ncol(m_y)
     K  <- length(l_x)

     v_psi   <- m_theta[, 1] # (N,)
     m_beta  <- m_theta[, 2:(K + 1), drop = FALSE] # (N,K), drop needed in case K == 1
     v_sgmsq <- m_theta[, K + 2] # (N,)

     v_sgm4h <- v_sgmsq^2
     v_sgm6h <- v_sgmsq^3

     # compute residuals
     m_beta_times_x <- matrix(0, N, TT) # 0s needed for recursive sum; (N,T) with generic element \vbeta_{i}'\vx_{it}
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          v_beta_k <- m_beta[, k] # (N,)
          m_beta_times_x_k <- m_x_k * v_beta_k #!! (N,T) "xHadamard" (N,1) = (N,T)
          m_beta_times_x <- m_beta_times_x + m_beta_times_x_k
     }
     m_psi_times_ys <- m_ys * v_psi #!! (N,T) "xHadamard" (N,1) = (N,T)
     m_eps <- m_y - m_psi_times_ys - m_beta_times_x # (N,T)
     v_ssr <- rowSums(m_eps^2) # (N,) sum of squared residuals: sum_t{(eps_it)^2} !!crossprod()??

     m_Psi <- diag(v_psi, nrow = N) # (N,N)
     m_A <- diag(N) - (m_Psi %*% m_W)

     det_mA <- det(m_A)
     if (det_mA <= 0) {
          stop("`I - Psi W` is not invertible. Check that the parameter space is defined correctly.")
     }

     # first derivative
     m_Q <- m_W %*% solve(m_A)

     # -------------------------------------------------------------------------- #
     # Inverse H matrix
     m_H11 <- (m_Q * t(m_Q)) + diag(rowSums(m_ys^2) / v_sgmsq / TT, nrow = N) # N x N
     m_H13 <- diag(rowSums(m_ys * m_eps) / v_sgm4h / TT) # N x N
     m_H33 <- diag(-(1 / 2 / v_sgm4h) + (v_ssr / v_sgm6h / TT)) # N x N

     m_H12  <- matrix(0, N, N * K)
     invH22 <- matrix(0, N * K, N * K)
     m_H23  <- matrix(0, N * K, N)

     for (i in 1:N) {
          ind <- ((i - 1) * K + 1):(i * K) #!! should I wrap this expression in c()?
          v_ysi <- m_ys[i, ] # (T,)
          m_Xi <- matrix(NA_real_, nrow = TT, ncol = K)
          for (k in 1:K) {
               m_x_k <- l_x[[k]] # (N,T)
               m_Xi[, k] <- m_x_k[i, ] # (T,)
          }
          v_epsi <- m_eps[i, ] # (T,)

          sgmsqi <- v_sgmsq[[i]]
          sgm4hi <- v_sgm4h[[i]]

          #stopifnot(K > 1)
          #!! 2021-09-21 I'm not 100% sure that this bit of code works as intended when K=1 
          m_H12[i, ind] <- (rbind(v_ysi) %*% m_Xi) / sgmsqi / TT # (1,K)
          invH22[ind, ind] <- solve(t(m_Xi) %*% m_Xi) * sgmsqi * TT # (K,K)
          m_H23[ind, i] <- (t(m_Xi) %*% cbind(v_epsi)) / sgm4hi / TT # (K,1)
     }

     m_Z11 <- m_H11
     m_Z12 <- cbind(m_H12, m_H13)

     invZ22 <- fn_inv_partitioned_a(invH22, m_H23, t(m_H23), m_H33)
     invH <- fn_inv_partitioned_b(m_Z11, m_Z12, t(m_Z12), invZ22)

     # -------------------------------------------------------------------------- #
     # J matrix (Equation~A.25)

     ## extract main diagonal
     v_q <- diag(m_Q) # (N,)
     m_q <- matrix(1, 1, TT) %x% cbind(v_q) # (N,T)

     m_sgmsq <- matrix(1, 1, TT) %x% cbind(v_sgmsq) # (N,T)
     m_sgm4h <- matrix(1, 1, TT) %x% cbind(v_sgm4h) # (N,T)
     m_dlogft_dvpsi <- (m_ys * m_eps / m_sgmsq) - m_q # (N,T)
     m_epssq <- m_eps^2
     m_dlogft_dvsgmsq <- (m_epssq / m_sgm4h / 2) - (1 / 2 / m_sgmsq) # (N,T)

     ##!! this part will probably be very slow because of the loop
     ## repmat: make v_sgmsq and m_eps as an array of dim (N,T,K) 
     ## allocate the same copy in the array, K times
     a_sgmsq <- array(NA_real_, c(N, TT, K))
     a_eps   <- array(NA_real_, c(N, TT, K))
     a_x     <- array(NA_real_, c(N, TT, K))
     for (k in 1:K) {
          a_eps  [, , k] <- m_eps # (N,T)
          a_sgmsq[, , k] <- m_sgmsq # (N,T)
          a_x    [, , k] <- l_x[[k]] # (N,T)
     }

     a_dlogft_dvbeta <- (a_eps * a_x) / a_sgmsq; # (N,T,K)
     a_dlogft_dvbeta_perm <- aperm(a_dlogft_dvbeta, c(3, 2, 1)); # (K,T,N)

     ## from (K,T,N) to (NK,T)
     m_dlogft_dvbeta <- matrix(NA_real_, K * N, TT); # (KN,T)
     for (i in 1:N) {
          v_ind <- seq.int(((i - 1) * K) + 1, (i * K)) # (K,)
          m_dlogft_dvbeta[v_ind, ] <- a_dlogft_dvbeta_perm[, , i]; # (K,T)
     }

     ## stack matrices on top of each others
     m_dlogft_dvtheta <- rbind(
          m_dlogft_dvpsi,
          m_dlogft_dvbeta,
          m_dlogft_dvsgmsq) # ((N+KN+N),T)

     ## J = DD' / T
     m_J <- (m_dlogft_dvtheta %*% t(m_dlogft_dvtheta)) / TT # (N(K+2),N(K+2))

     # -------------------------------------------------------------------------- #
     # standard variance
     v_var <- diag(invH) / TT # (N(K+2),)
     m_standard <- fn_reshape_theta_vec2mat(v_var, N, K) # (N,K+2)

     # -------------------------------------------------------------------------- #
     # sandwich variance
     m_invH_J_invH <- invH %*% m_J %*% invH # (N(K+2),N(K+2))

     ## extract the main diagonal
     v_var <- diag(m_invH_J_invH) / TT # (N(K+2),)

     ## reshape
     m_sandwich <- fn_reshape_theta_vec2mat(v_var, N, K) # (N,K+2)

     # -------------------------------------------------------------------------- #
     # return
     list(standard = m_standard, sandwich = m_sandwich)
}

# -------------------------------------------------------------------------- #
# reshape a vector of length N(K+2) into a matrix of order (N,K+2) 
fn_reshape_theta_vec2mat <- function(vec, N, K) {
     mat <- matrix(NA_real_, N, K + 2)
     mat[, 1] <- vec[1:N]
     mat[, 2:(K + 1)] <- t(matrix(vec[(N + 1):(N + (K * N))], K, N)) # (N,K)
     mat[, K + 2] <- vec[(N + (K * N) + 1):(N + (K * N) + N)] # (N,)

     # return
     mat
}

# -------------------------------------------------------------------------- #
fn_inv_partitioned_a <- function(invA, m_B, m_C, m_D) {
     m_C_invA <- m_C %*% invA
     m_E <- m_D - (m_C_invA %*% m_B)
     invE <- solve(m_E)
     m_invA_B_invE <- invA %*% m_B %*% invE

     invH11 <- invA + (m_invA_B_invE %*% m_C_invA)
     invH12 <- -m_invA_B_invE
     invH21 <- -invE %*% m_C_invA
     invH22 <- invE

     invH <- rbind(
          cbind(invH11, invH12),
          cbind(invH21, invH22)
          )

     # return
     invH
}

# -------------------------------------------------------------------------- #
fn_inv_partitioned_b <- function(m_A, m_B, m_C, invD) {
     m_B_invD <- m_B %*% invD
     m_F <- m_A - (m_B_invD %*% m_C)
     invF <- solve(m_F)
     m_invD_C_invF <- invD %*% m_C %*% invF

     invH11 <- invF
     invH12 <- -invF %*% m_B_invD
     invH21 <- -m_invD_C_invF
     invH22 <- invD + (m_invD_C_invF %*% m_B_invD)

     invH <- rbind(
          cbind(invH11, invH12),
          cbind(invH21, invH22)
          )

     # return
     invH
}
