% Predicting Overall Survival Using Early PFS Evidence
% Jinjie Chen
% Internal Methodology Seminar

# 1. Motivation

- Overall Survival (OS) is the gold standard endpoint in oncology.
- However OS often matures much later than PFS.

Typical timeline:

- PFS readout: early
- OS readout: delayed

Goal:

Estimate the probability that the final OS result will be successful.

---

# 2. Decision Quantity

Define success as an OS hazard ratio threshold.

$$
HR_{target}
$$

On the log scale

$$
\ell_{target} = \log(HR_{target})
$$

Predictive probability of success:

$$
PoS = P(\theta_{OS}^{(c)} < \ell_{target})
$$

---

# 3. Available Data

Two types of information:

Historical trials

Current trial

Historical trials provide the relationship between OS and PFS.

Current trial provides early evidence.

---

# 4. Historical Trial Data

For each historical study we observe

OS effect estimate

$$
\hat{\theta}_{OS,i}
$$

PFS effect estimate

$$
\hat{\theta}_{PFS,i}
$$

Standard errors

$$
SE_{OS,i}, \quad SE_{PFS,i}
$$

Within study correlation

$$
\rho_i
$$

---

# 5. Vector Representation

Define observed estimates

$$
y_i =
\begin{pmatrix}
\hat{\theta}_{OS,i} \\
\hat{\theta}_{PFS,i}
\end{pmatrix}
$$

True treatment effects

$$
\theta_i =
\begin{pmatrix}
\theta_{OS,i} \\
\theta_{PFS,i}
\end{pmatrix}
$$

---

# 6. Within Study Sampling Model

Assume asymptotic normality

$$
y_i | \theta_i \sim N(\theta_i, S_i)
$$

Within study covariance matrix

$$
S_i =
\begin{pmatrix}
SE_{OS,i}^2 & \rho_i SE_{OS,i} SE_{PFS,i} \\
\rho_i SE_{OS,i} SE_{PFS,i} & SE_{PFS,i}^2
\end{pmatrix}
$$

---

# 7. Between Study Model

Assume random effects

$$
\theta_i \sim N(\eta, \Sigma_0)
$$

where

$$
\eta =
\begin{pmatrix}
\eta_{OS} \\
\eta_{PFS}
\end{pmatrix}
$$

---

# 8. Between Study Covariance

$$
\Sigma_0 =
\begin{pmatrix}
\Sigma_{OO} & \Sigma_{OP} \\
\Sigma_{OP} & \Sigma_{PP}
\end{pmatrix}
$$

This matrix captures

- heterogeneity across trials
- correlation between OS and PFS effects

---

# 9. Marginal Distribution

Combining levels

$$
y_i \sim N(\eta, S_i + \Sigma_0)
$$

We estimate

$$
\eta , \Sigma_0
$$

using multivariate meta analysis.

---

# 10. Estimation

Two approaches

REML

Method of Moments

Primary analysis uses REML.

Method of Moments is used as sensitivity analysis.

---

# 11. Current Trial Evidence

Current trial provides

PFS estimate

$$
\hat{\theta}_{PFS}^{(c)}
$$

Standard error

$$
SE_{PFS}^{(c)}
$$

Goal:

Predict the OS treatment effect.

---

# 12. Conditional Prediction

Using multivariate normal theory

$$
E(\theta_{OS}^{(c)}|\theta_{PFS}^{(c)}) =
\eta_{OS} +
\frac{\Sigma_{OP}}{\Sigma_{PP}}
(\theta_{PFS}^{(c)}-\eta_{PFS})
$$

---

# 13. Conditional Variance

$$
Var(\theta_{OS}^{(c)}|\theta_{PFS}^{(c)}) =
\Sigma_{OO} -
\frac{\Sigma_{OP}^2}{\Sigma_{PP}}
$$

---

# 14. Predictive Distribution

Accounting for estimation uncertainty

$$
\theta_{OS}^{(c)} \sim N(\mu_{pred}, \sigma_{pred}^2)
$$

---

# 15. Probability of Success

Closed form probability

$$
PoS =
\Phi
\left(
\frac{\ell_{target}-\mu_{pred}}
{\sigma_{pred}}
\right)
$$

---

# 16. OS Interim Information

Sometimes interim OS is available.

Observed interim estimate

$$
\hat{\theta}_{OS}^{int}
$$

Standard error

$$
SE_{OS}^{int}
$$

---

# 17. Joint Updating

Observed data

$$
y =
\begin{pmatrix}
\hat{\theta}_{OS}^{int} \\
\hat{\theta}_{PFS}^{(c)}
\end{pmatrix}
$$

---

# 18. Posterior Distribution

Prior

$$
\theta^{(c)} \sim N(\eta,\Sigma_B)
$$

Posterior

$$
\theta^{(c)}|y \sim N(m,V)
$$

---

# 19. Posterior Parameters

$$
V = (\Sigma_B^{-1}+S^{-1})^{-1}
$$

$$
m = V(\Sigma_B^{-1}\eta + S^{-1}y)
$$

---

# 20. Summary

Framework

1 Historical evidence informs OS–PFS relationship

2 Current PFS predicts OS

3 Interim OS can refine prediction

4 PoS provides a quantitative decision metric