## Version 1.0.0
First release

## Changes in Version 1.0.1
- As shown in the paper, we propose a simulated general method of moments (SGMM) for the SAR model (see function smmSAR and Section 2 of our vignette).
- We can now estimate the maximal bias of the instrumental variable estimator (see Section 1.1 and 1.2 of our vignette).

## Changes in Version 1.0.2
- We provide a smoother simulator of adjacency matrices in the SGMM approach.
- We add weights to the probit/logit network formation model.
- We allows the use of an initial probit/logit estimate of  $\rho$, where the observed part of the network is assumed non-stochastic in the MCMC. This is a quite different from using an initial probit/logit estimate as prior distribution of $\rho$. In this latter case, $\rho$ is updated using, among others, information from the observed part of the network. In the first case, $\rho$ and the unobserved part of the network are updated using information in $y$, where the initial estimate acts as prior distribution of $\rho$. Information from the observed part of the network is not used to update $\rho$. This information is included in the initial estimate.

## Changes in Version 1.0.3
Adjustments with Eigen 3.4

## Changes in Version 1.0.4
Adjustments with CDatanet 2.2.0

## Changes in Version 1.1.0
- `remove.ids` function is added.
- The moment function in `smmSAR` when `GX` is observed is corrected following the revision of the paper in March 2025.

## Changes in Version 1.1.1
- Accommodate changes in Armadillo 13