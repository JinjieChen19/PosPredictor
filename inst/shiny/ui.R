library(shiny)

navbarPage(
  title = "PosPredictor: Using PFS to Predict OS",
  theme = NULL,
  id    = "main_nav",

  # ===========================================================================
  # Tab 1: Model Description
  # ===========================================================================
  tabPanel(
    "Model Description",
    fluidPage(
      h2("Hierarchical Random-Effects Meta-Analysis Model"),
      p("This application computes the posterior Probability of Success (PoS)",
        "for an ongoing OS readout, borrowing strength from historical trials",
        "through a three-level hierarchical (random-effects) meta-analytic model."),

      hr(),
      h3("Model Structure"),

      h4("Level 3 — Population Hyperparameters"),
      withMathJax(
        p("The population mean log(HR) vector follows a bivariate normal prior:"),
        p("$$\\boldsymbol{\\mu} = (\\mu_{OS}, \\mu_{PFS})^\\top, \\quad
           \\mu_{OS} \\sim N(m_{OS}, s_{OS}^2), \\quad
           \\mu_{PFS} \\sim N(m_{PFS}, s_{PFS}^2)$$"),
        p("The between-trial covariance matrix is:"),
        p("$$\\boldsymbol{\\Sigma} = \\begin{pmatrix}
             \\tau_{OS}^2 & \\rho\\,\\tau_{OS}\\tau_{PFS} \\\\
             \\rho\\,\\tau_{OS}\\tau_{PFS} & \\tau_{PFS}^2
           \\end{pmatrix}$$"),
        p("Priors: \\(\\tau_{OS}, \\tau_{PFS} \\sim \\text{half-Normal}\\);",
          " \\(\\rho\\) is modelled via Fisher-z transform: \\(\\tanh^{-1}(\\rho) \\sim N(0, \\sigma_z^2)\\).")
      ),

      h4("Level 2 — Trial-Specific Effects"),
      withMathJax(
        p("For each trial \\(k = 1, \\ldots, K\\) (historical) and the current trial:"),
        p("$$\\boldsymbol{\\theta}_k = (\\theta_{k,OS}, \\theta_{k,PFS})^\\top
           \\sim \\text{MVN}(\\boldsymbol{\\mu}, \\boldsymbol{\\Sigma})$$"),
        p("These represent the ", strong("true"), " treatment effects in each trial."),
        p(strong("Non-Centered Parameterisation"), "(for efficient MCMC):"),
        p("$$\\boldsymbol{\\theta}_{raw,k} \\sim N(\\mathbf{0}, I), \\quad
           \\boldsymbol{\\theta}_k = \\boldsymbol{\\mu} + L\\,\\boldsymbol{\\theta}_{raw,k}$$"),
        p("where \\(L\\) is the Cholesky factor of \\(\\boldsymbol{\\Sigma}\\).")
      ),

      h4("Level 1 — Observed Data (Likelihood)"),
      withMathJax(
        p("Observed log(HR) pairs \\(\\mathbf{y}_k = (y_{k,OS}, y_{k,PFS})^\\top\\) are",
          " noisy measurements of the true effects:"),
        p("$$\\mathbf{y}_k \\sim \\text{MVN}(\\boldsymbol{\\theta}_k, W_k)$$"),
        p("where \\(W_k\\) is the known within-trial covariance matrix derived",
          " from observed standard errors and within-trial correlation.")
      ),

      h4("Probability of Success"),
      withMathJax(
        p("$$\\text{PoS}_{OS} = \\Pr(\\theta_{current,OS} < \\text{target}_{OS} \\mid \\text{data})$$"),
        p("Example: if \\(\\text{target}_{OS} = \\log(0.74) \\approx -0.30\\) and",
          " PoS = 0.67, there is a 67% probability the final trial will demonstrate",
          " at least a 26% risk reduction in OS.")
      ),

      hr(),
      h3("Prior Specifications"),
      tableOutput("prior_table"),

      hr(),
      h3("Information Flow"),
      tags$ul(
        tags$li("Population parameters (μ, Σ) govern the distribution of trial-specific effects."),
        tags$li("Trial-specific effects (θ_k) represent the true treatment effects in each trial."),
        tags$li("Observed data (y_k) provide noisy measurements of the true effects."),
        tags$li("The model 'borrows strength' across trials while accounting for heterogeneity.")
      ),

      hr(),
      h3("Computational Approaches"),
      tags$ul(
        tags$li(strong("Closed-Form (Tab 3):"),
                " Gaussian conditioning on the observed PFS log-HR.",
                " The empirical joint covariance of (OS, PFS) from historical trials is",
                " used to derive the posterior mean and variance of OS given the current",
                " PFS readout. No between-trial heterogeneity parameters need to be specified."),
        tags$li(strong("Full Bayesian MCMC (Tab 4):"),
                " Hamiltonian Monte Carlo via Stan (2000 iterations, 1000 warmup,",
                " 4 chains, adapt_delta = 0.97). Stan model is",
                strong("pre-compiled"), "and saved as .RDS for instant reload.",
                " All priors passed as data — no recompilation needed.")
      )
    )
  ),

  # ===========================================================================
  # Tab 2: Data Input
  # ===========================================================================
  tabPanel(
    "Data Input",
    fluidPage(
      h2("Trial Data Input"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Historical Data Simulation"),
          sliderInput("K", "Number of historical trials", 5, 40, 27, step = 1),
          sliderInput("mu_os_true", "True μ_OS (pop. mean log HR OS)",
                      -0.8, 0.2, -0.30, step = 0.01),
          sliderInput("mu_pfs_true", "True μ_PFS (pop. mean log HR PFS)",
                      -0.8, 0.2, -0.45, step = 0.01),
          sliderInput("tau_os_true", "True τ_OS (between-trial SD)",
                      0.01, 0.5, 0.15, step = 0.01),
          sliderInput("tau_pfs_true", "True τ_PFS (between-trial SD)",
                      0.01, 0.5, 0.12, step = 0.01),
          sliderInput("rho_true", "True ρ (between-trial correlation)",
                      -0.95, 0.95, 0.65, step = 0.05),
          numericInput("sim_seed", "Random seed", 20260212, min = 1),
          actionButton("simulate_btn", "Simulate Historical Data",
                       class = "btn-primary"),
          hr(),
          h4("Current Trial"),
          numericInput("cur_y_os",  "Current log(HR) OS",  -0.30, step = 0.01),
          numericInput("cur_y_pfs", "Current log(HR) PFS", -0.60, step = 0.01),
          numericInput("cur_se_os",  "Current SE (OS)",   0.20, min = 0.01, step = 0.01),
          numericInput("cur_se_pfs", "Current SE (PFS)",  0.10, min = 0.01, step = 0.01),
          sliderInput("cur_within_corr", "Current within-trial corr.",
                      0.2, 0.95, 0.65, step = 0.05),
          numericInput("target_os", "Success threshold (target log HR OS)",
                       -0.30, step = 0.01)
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel(
              "Scatter Plot",
              br(),
              plotOutput("scatter_plot", height = "500px"),
              br(),
              p("Circle size is proportional to number of OS events.",
                "The regression line (weighted OLS) is shown in red.",
                "The dashed line represents perfect surrogacy (slope = 1).",
                "The orange cross marks the current trial interim readout.")
            ),
            tabPanel(
              "Data Table",
              br(),
              DT::dataTableOutput("hist_table")
            )
          )
        )
      )
    )
  ),

  # ===========================================================================
  # Tab 3: Closed-Form PoS
  # ===========================================================================
  tabPanel(
    "Closed-Form PoS",
    fluidPage(
      h2("Closed-Form PoS via Gaussian Conditioning on PFS"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Population Mean Assumptions"),
          p(style = "font-size:12px; color:#555;",
            "These are the assumed population means for OS and PFS log(HR).",
            "The empirical covariance is estimated directly from the historical data."),
          numericInput("cf_mu_os_prior_mean",  "Population mean μ_OS",  -0.30, step = 0.01),
          numericInput("cf_mu_pfs_prior_mean", "Population mean μ_PFS", -0.45, step = 0.01),
          actionButton("run_cf_btn", "Compute Closed-Form PoS",
                       class = "btn-success")
        ),
        mainPanel(
          width = 9,
          uiOutput("cf_pos_box"),
          br(),
          plotOutput("cf_posterior_os_plot", height = "350px"),
          br(),
          h4("Formula Summary"),
          withMathJax(
            p("Using empirical sample covariance \\(\\hat{S}\\) from historical log-HR pairs:"),
            p("$$\\hat{m} = \\mu_{OS} + \\frac{\\hat{S}_{12}}{\\hat{S}_{22} + \\sigma_{PFS}^2}(y_{PFS} - \\mu_{PFS})$$"),
            p("$$\\hat{v} = \\hat{S}_{11} - \\frac{\\hat{S}_{12}^2}{\\hat{S}_{22} + \\sigma_{PFS}^2}$$"),
            p("$$\\text{PoS} = \\Phi\\!\\left(\\frac{\\text{target}_{OS} - \\hat{m}}{\\sqrt{\\hat{v} + \\sigma_{PFS}^2}}\\right)$$"),
            p(style = "font-size:11px; color:#777;",
              "Note: following the reference simulation, σ_PFS (SE of current PFS readout)",
              "is used in the z-score denominator as the residual prediction uncertainty term.")
          )
        )
      )
    )
  ),

  # ===========================================================================
  # Tab 4: MCMC PoS
  # ===========================================================================
  tabPanel(
    "MCMC PoS (Stan)",
    fluidPage(
      h2("Full Bayesian MCMC Posterior PoS via Stan"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Prior Hyperparameters"),
          numericInput("mc_mu_os_mean",  "Prior mean μ_OS",  -0.30, step = 0.01),
          numericInput("mc_mu_os_sd",    "Prior SD μ_OS",     0.50, min = 0.01, step = 0.05),
          numericInput("mc_mu_pfs_mean", "Prior mean μ_PFS", -0.45, step = 0.01),
          numericInput("mc_mu_pfs_sd",   "Prior SD μ_PFS",    0.50, min = 0.01, step = 0.05),
          numericInput("mc_tau_os_sd",   "τ_OS half-Normal SD",  0.25, min = 0.01, step = 0.05),
          numericInput("mc_tau_pfs_sd",  "τ_PFS half-Normal SD", 0.25, min = 0.01, step = 0.05),
          numericInput("mc_rho_z_sd",    "Fisher-z(ρ) prior SD", 1.50, min = 0.10, step = 0.10),
          h4("MCMC Settings"),
          numericInput("mc_iter",    "Total iterations",  2000, min = 500,  step = 500),
          numericInput("mc_warmup",  "Warmup iterations", 1000, min = 250,  step = 250),
          numericInput("mc_chains",  "Chains",               4, min = 1,    step = 1),
          numericInput("mc_adapt",   "Adapt delta",       0.97, min = 0.80, max = 0.999, step = 0.01),
          numericInput("mc_seed",    "Random seed",        123, min = 1),
          h4("Model Compilation"),
          p(style = "font-size:12px; color:#555;",
            "The compiled model is cached automatically in the R user cache directory and reused across sessions without recompilation."),
          actionButton("compile_btn", "Compile & Save Stan Model",
                       class = "btn-warning"),
          br(), br(),
          actionButton("run_mcmc_btn", "Run MCMC",
                       class = "btn-danger")
        ),
        mainPanel(
          width = 9,
          uiOutput("mcmc_pos_box"),
          br(),
          uiOutput("mcmc_convergence_ui"),
          br(),
          tabsetPanel(
            tabPanel(
              "Trace Plots",
              br(),
              plotOutput("mcmc_trace_plot", height = "500px")
            ),
            tabPanel(
              "Posterior Densities",
              br(),
              plotOutput("mcmc_density_plot", height = "500px")
            ),
            tabPanel(
              "Pairs Plot",
              br(),
              plotOutput("mcmc_pairs_plot", height = "500px")
            ),
            tabPanel(
              "Diagnostics Table",
              br(),
              tableOutput("mcmc_diag_table")
            )
          )
        )
      )
    )
  )
)
