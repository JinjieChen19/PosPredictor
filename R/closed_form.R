#' Compute PoS via Closed-Form Gaussian Conditioning on PFS
#'
#' Uses the empirical joint distribution of historical OS and PFS log-HRs
#' together with a Gaussian conditioning formula to compute the posterior
#' probability of success (PoS) for the current trial's OS endpoint, given
#' the observed PFS log-HR.
#'
#' Algorithm (mirrors the reference simulation approach):
#' \enumerate{
#'   \item Compute the sample covariance matrix \eqn{S} of the historical
#'         observed log-HR pairs \eqn{(y_{OS}, y_{PFS})}.
#'   \item Bayesian linear prediction:
#'         \eqn{\hat{m} = \mu_{OS} + S_{12}/(S_{22} + \sigma_{PFS}^2) \times (y_{PFS,cur} - \mu_{PFS})}
#'   \item Conditional variance:
#'         \eqn{\hat{v} = S_{11} - S_{12}^2 / (S_{22} + \sigma_{PFS}^2)}
#'   \item PoS: \eqn{\Phi\!\left(\frac{\text{target}_{OS} - \hat{m}}{\sqrt{\hat{v} + \sigma_{PFS}^2}}\right)}
#' }
#'
#' \strong{Note on z-score denominator}: following the reference simulation
#' code, the effective SD uses \eqn{\sigma_{PFS}^2} (not \eqn{\sigma_{OS}^2})
#' to account for the residual prediction uncertainty from the PFS measurement.
#' \code{current_se_os} is accepted but not used in the formula; it is retained
#' in the signature for forward-compatibility.
#'
#' @param hist_data data.frame from \code{\link{simulate_historical_data}}.
#' @param current_y_pfs Numeric. Current trial observed log(HR) for PFS.
#' @param current_se_pfs Numeric. Current trial SE for PFS log(HR).
#' @param target_os Numeric. Success threshold for OS log(HR) (e.g. -0.30).
#' @param mu_os Numeric. Assumed population mean log(HR) for OS (default -0.30).
#' @param mu_pfs Numeric. Assumed population mean log(HR) for PFS (default -0.45).
#' @param current_se_os Numeric. Current trial SE for OS (retained for
#'   forward-compatibility; not used in the current formula).
#' @return A named list with: pos, m_post, v_post, sd_eff,
#'   credible_interval_os.
#' @export
compute_pos_closed_form <- function(hist_data,
                                     current_y_pfs,
                                     current_se_pfs,
                                     target_os     = -0.30,
                                     mu_os         = -0.30,
                                     mu_pfs        = -0.45,
                                     current_se_os = NULL) {

  # Sample covariance from historical observed log-HR pairs
  S_obs      <- stats::cov(data.frame(OS  = hist_data$log_hr_os,
                                      PFS = hist_data$log_hr_pfs))
  var_pfs    <- S_obs[2, 2]
  var_os     <- S_obs[1, 1]
  cov_os_pfs <- S_obs[1, 2]

  # Gaussian conditioning: E[theta_OS | y_PFS] and Var(theta_OS | y_PFS)
  denom  <- var_pfs + current_se_pfs^2
  m_post <- mu_os + cov_os_pfs / denom * (current_y_pfs - mu_pfs)
  v_post <- max(var_os - cov_os_pfs^2 / denom, 0)

  # Effective SD used for z-score (matches reference formula)
  sd_eff <- sqrt(v_post + current_se_pfs^2)

  # PoS = P(theta_OS < target_os)
  z_score <- (target_os - m_post) / sd_eff
  pos     <- stats::pnorm(z_score)

  # 95% credible interval for theta_OS (on the same predictive distribution)
  ci_os <- stats::qnorm(c(0.025, 0.975), mean = m_post, sd = sd_eff)

  list(
    pos                  = pos,
    m_post               = m_post,
    v_post               = v_post,
    sd_eff               = sd_eff,
    credible_interval_os = ci_os
  )
}
