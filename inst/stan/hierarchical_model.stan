// Hierarchical Random-Effects Meta-Analysis for PFS → OS Surrogacy
//
// Three-Level Model:
//   Level 3: Population parameters  mu ~ MVN(mu_prior, diag(sigma_prior^2))
//            tau_os, tau_pfs ~ half-Normal;  rho via Fisher-z transform
//   Level 2: Trial-specific effects theta_k ~ MVN(mu, Sigma)
//            Non-centered parameterisation for efficient MCMC
//   Level 1: Likelihood  y_k ~ MVN(theta_k, W_k) where y_k are observed
//            log(HR) pairs and W_k is the known within-trial covariance
//
// Priors are passed as data so the model does NOT need recompilation
// when prior hyperparameters or computation settings change.

data {
  int<lower=1> K;                        // number of historical trials

  // Historical trial data (2 x K: row 1 = OS, row 2 = PFS)
  array[K] vector[2] y_hist;            // observed log HRs
  array[K] matrix[2, 2] W_hist;        // within-trial covariance matrices

  // Current trial
  vector[2] y_current;                   // observed log HRs (OS, PFS)
  matrix[2, 2] W_current;               // within-trial covariance

  real target_os;                        // OS success threshold (e.g. log(0.74))

  // Prior hyperparameters (passed as data – no recompilation needed)
  real    mu_os_prior_mean;
  real<lower=0> mu_os_prior_sd;
  real    mu_pfs_prior_mean;
  real<lower=0> mu_pfs_prior_sd;
  real<lower=0> tau_os_prior_sd;        // half-Normal SD for tau_os
  real<lower=0> tau_pfs_prior_sd;       // half-Normal SD for tau_pfs
  // rho_z is the unconstrained Fisher transform of the between-trial correlation.
  // rho = tanh(rho_z) maps rho_z back to (-1, 1).
  // Setting rho_z_prior_mean = atanh(0.65) ≈ 0.775 centers the prior around
  // correlation 0.65; rho_z_prior_mean = 0 gives a symmetric prior around rho = 0.
  real rho_z_prior_mean;                 // prior mean for Fisher-z(rho)
  real<lower=0> rho_z_prior_sd;         // prior SD  for Fisher-z(rho)
}

parameters {
  // Level 3: Population parameters
  vector[2] mu;                          // (mu_os, mu_pfs)
  real<lower=0> tau_os;                 // between-trial SD for OS
  real<lower=0> tau_pfs;               // between-trial SD for PFS
  real rho_z;                            // Fisher-z transform of rho

  // Level 2: Non-centered trial-specific effects
  array[K] vector[2] theta_raw_hist;    // standard normal auxiliaries
  vector[2] theta_raw_current;           // for current trial
}

transformed parameters {
  real<lower=-1, upper=1> rho;
  matrix[2, 2] Sigma;
  matrix[2, 2] L_Sigma;
  array[K] vector[2] theta_hist;
  vector[2] theta_current;

  // Between-trial covariance and its Cholesky factor
  rho         = tanh(rho_z);
  Sigma[1, 1] = square(tau_os);
  Sigma[2, 2] = square(tau_pfs);
  Sigma[1, 2] = rho * tau_os * tau_pfs;
  Sigma[2, 1] = Sigma[1, 2];
  L_Sigma     = cholesky_decompose(Sigma);

  // Non-centered parameterisation: theta_k = mu + L * theta_raw_k
  for (k in 1:K) {
    theta_hist[k] = mu + L_Sigma * theta_raw_hist[k];
  }
  theta_current = mu + L_Sigma * theta_raw_current;
}

model {
  // --- Level 3 priors ---
  mu[1]   ~ normal(mu_os_prior_mean,  mu_os_prior_sd);
  mu[2]   ~ normal(mu_pfs_prior_mean, mu_pfs_prior_sd);
  tau_os  ~ normal(0, tau_os_prior_sd);     // half-normal (tau_os >= 0)
  tau_pfs ~ normal(0, tau_pfs_prior_sd);    // half-normal
  // rho_z ~ N(rho_z_prior_mean, rho_z_prior_sd): Fisher-z prior for correlation.
  // rho = tanh(rho_z) constrains the between-trial correlation to (-1, 1).
  rho_z   ~ normal(rho_z_prior_mean, rho_z_prior_sd);

  // --- Level 2: non-centered standard-normal auxiliaries ---
  for (k in 1:K) {
    theta_raw_hist[k] ~ std_normal();
  }
  theta_raw_current ~ std_normal();

  // --- Level 1: likelihood ---
  for (k in 1:K) {
    y_hist[k] ~ multi_normal(theta_hist[k], W_hist[k]);
  }
  y_current ~ multi_normal(theta_current, W_current);
}

generated quantities {
  // PoS = Pr(theta_current_OS < target_os | data)
  // step(x) = 1 if x >= 0, so step(target_os - theta_current[1]) = 1 iff
  // theta_current[1] <= target_os, i.e. the OS effect meets the threshold.
  real pos_os = step(target_os - theta_current[1]);

  // Diagnostic aliases: copies of transformed parameters for user-facing
  // posterior inspection (e.g. in bayesplot).  They do not change inference.
  real rho_out          = rho;               // between-trial correlation on (-1, 1)
  real theta_os_current  = theta_current[1]; // current trial true OS log(HR)
  real theta_pfs_current = theta_current[2]; // current trial true PFS log(HR)
}
