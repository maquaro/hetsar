#' @export
summary.hetsar <- function(l_results, MG = FALSE, compact = TRUE, robust = TRUE) {

     if (!inherits(l_results, "hetsar")) {
          stop("input needs to be of class 'hetsar'")
     }

     if (!MG) {
          m_theta       <- l_results[["coefficients"]]
          m_se_sandwich <- l_results[["Std.err.sandwich"]]
          m_se_standard <- l_results[["Std.err.standard"]]
          T. <- l_results[["T"]]
          N. <- nrow(m_theta)

          if (robust) {
               m_se <- m_se_sandwich
          } else {
               m_se <- m_se_standard
          }

          m_tval  <- m_theta / m_se
          m_atval <- abs(m_tval)
          m_stars <- (m_atval > 1.645) + (m_atval > 1.96) + (m_atval > 2.33)
          m_95l <- m_theta - (1.96 * m_se)
          m_95u <- m_theta + (1.96 * m_se)

          v_id <- rownames(m_theta)
          v_name <- colnames(m_theta)
          K <- length(v_name)

          if (N. < 7) { #cc#
               compact <- TRUE
               i1 <- Inf
               i2 <- -Inf
          } else {
               i1 <- 4
               i2 <- N. - 3
          }
          cat("SAR with heterogeneous coefficients\n\n")
          cat("Number of obs = ")
          cat(format(N. * T., big.mark = ",", scientific = FALSE))
          cat(sprintf(" (N = %d, T = %d)\n", N., T.))
          cat(sprintf('%s\n', paste(rep('=', 78), collapse = "")))
          cat(sprintf('%15s', 'id'))
          cat(sprintf(' %13s', 'Coef.'))
          cat(sprintf(' %12s', 'Std. Err.'))
          cat(sprintf(' %8s', 'z'))
          cat(sprintf(' %25s\n', '[95% Conf. Interval]'))
          cat(sprintf('%s\n', paste(rep('-', 78), collapse = "")))
          for (k in 1:K) {
               cat(sprintf(' %s\n', v_name[k]))
               for (i in 1:N.) {
                    if (i == i1) {
                         cat(sprintf('% 17s |', "\u22ee"))
                         cat(sprintf(' % 10s',  "\u22ee"))
                         cat(sprintf('%-3s', ' '))
                         cat(sprintf(' %12s',   "\u22ee"))
                         cat(sprintf(' %12s',   "\u22ee"))
                         cat(sprintf(' %12s',   "\u22ee"))
                         cat(sprintf(' %12s\n', "\u22ee"))
                    } else if (i < i1 || i > i2) {
                         id = v_id[i]
                         if (typeof(id) == "character") {
                              cat(sprintf('% 15s |', id)) #!! TODO what it names are too long?
                         } else {
                              cat(sprintf('% 15d |', id))
                         }

                         cat(sprintf(' % 10.3f', m_theta[i, k])) #!! TODO what if these numbers are too small / too large? How is done in Stata
                         if        (m_stars[i, k] == 0) { cat(sprintf('%-3s', ' '))
                         } else if (m_stars[i, k] == 1) { cat(sprintf('%-3s', '*'))
                         } else if (m_stars[i, k] == 2) { cat(sprintf('%-3s', '**'))
                         } else if (m_stars[i, k] == 3) { cat(sprintf('%-3s', '***'))
                         }
                         cat(sprintf(' %10.3f', m_se[i, k]))
                         cat(sprintf(' % 10.3f', m_tval[i, k]))
                         cat(sprintf(' % 10.3f',   m_95l[i, k]))
                         cat(sprintf(' % 10.3f\n', m_95u[i, k]))
                    }
               }
               if (k < K) {
                    cat(sprintf('%s\n', paste(rep('-', 78), collapse = "")))
               } else {
                    cat(sprintf('%s\n', paste(rep('=', 78), collapse = "")))
               }
          }

          # pack all results into a data-frame
          df_hetsar <- data.frame(
               variable = rep(v_name, each = N.),
               id = rep(v_id, times = K),
               coef = c(m_theta),
               stars = c(m_stars),
               se = c(m_se),
               tval = c(m_tval),
               ci95l = c(m_95l),
               ci95u = c(m_95u)
          )

          return(invisible(df_hetsar))

     } else {
          m_theta <- l_results[["coefficients"]]
          T. <- l_results[["T"]]
          N. <- nrow(m_theta)
          K <- ncol(m_theta)

          # remove sgmsq
          v_name <- colnames(m_theta)[1:(K - 1)]
          K <- K - 1

          df_mg <- data.frame(
               variable = v_name,
               coef  = rep(NA_real_, K),
               stars = rep(NA_real_, K),
               se    = rep(NA_real_, K),
               tval  = rep(NA_real_, K),
               ci95l = rep(NA_real_, K),
               ci95u = rep(NA_real_, K))

          for (k in 1:K) {
               name_k <- v_name[k]
               vec <- m_theta[, name_k] # (N,)
               l_out <- fn_mg(vec)
               df_mg[k, 'coef' ] <- l_out$MG_coef
               df_mg[k, 'stars'] <- l_out$MG_stars
               df_mg[k, 'se'   ] <- l_out$MG_se
               df_mg[k, 'tval' ] <- l_out$MG_tval
               df_mg[k, 'ci95l'] <- l_out$MG_95l
               df_mg[k, 'ci95u'] <- l_out$MG_95u
          }

          cat("SAR with heterogeneous coefficients\n")
          cat("Mean-group estimates\n\n")
          cat("Number of obs = ")
          cat(format(N. * T., big.mark = ",", scientific = FALSE))
          cat(sprintf(" (N = %d, T = %d)\n", N., T.))
          cat(sprintf('%s\n', paste(rep('=', 78), collapse = "")))
          cat(sprintf('%15s', 'variable'))
          cat(sprintf(' %13s', 'Coef.'))
          cat(sprintf(' %12s', 'Std. Err.'))
          cat(sprintf(' %8s', 'z'))
          cat(sprintf(' %25s\n', '[95% Conf. Interval]'))
          cat(sprintf('%s+%s\n', paste(rep('-', 16), collapse = ""), paste(rep('-', 61), collapse = "")))
          for (k in 1:K) {
               cat(sprintf('% 15s |', df_mg[k, 'variable']))

               cat(sprintf(' % 10.3f', df_mg[k, 'coef']))
               if        (df_mg[k, 'stars'] == 0) { cat(sprintf('%-3s', ' '))
               } else if (df_mg[k, 'stars'] == 1) { cat(sprintf('%-3s', '*'))
               } else if (df_mg[k, 'stars'] == 2) { cat(sprintf('%-3s', '**'))
               } else if (df_mg[k, 'stars'] == 3) { cat(sprintf('%-3s', '***'))
               }
               cat(sprintf(' %10.3f',    df_mg[k, 'se'   ]))
               cat(sprintf(' % 10.3f',   df_mg[k, 'tval' ]))
               cat(sprintf(' % 10.3f',   df_mg[k, 'ci95l']))
               cat(sprintf(' % 10.3f\n', df_mg[k, 'ci95u']))
          }
          if (k < K) {
               cat(sprintf('%s\n', paste(rep('-', 78), collapse = "")))
          } else {
               cat(sprintf('%s\n', paste(rep('=', 78), collapse = "")))
          }

          return(invisible(df_mg))
     }
}

fn_mg <- function(v_psi) {
     N <- length(v_psi)
     psiMG <- mean(v_psi)
     v_psi_exc <- v_psi - psiMG # (N,1)
     v_psi_excsq <- v_psi_exc * v_psi_exc
     somma <- sum(v_psi_excsq)
     varianza <- somma / N / (N - 1)
     psiMG_se <- sqrt(varianza)
     psiMG_tval <- psiMG / psiMG_se
     psiMG_atval <- abs(psiMG_tval)
     psiMG_stars <- (psiMG_atval > 1.645) + (psiMG_atval > 1.96) + (psiMG_atval > 2.33)
     psiMG_95l <- psiMG - (1.96 * psiMG_se)
     psiMG_95u <- psiMG + (1.96 * psiMG_se)

     # save these three scalars into a struc
     l_return          <- list()
     l_return$MG_coef  <- psiMG
     l_return$MG_stars <- psiMG_stars
     l_return$MG_se    <- psiMG_se
     l_return$MG_tval  <- psiMG_tval
     l_return$MG_95l   <- psiMG_95l
     l_return$MG_95u   <- psiMG_95u

     l_return
}
