------------------------------------------- CPP FUNCTIONS --------------------------------------------
File CppFunctions.cpp contains C++ source functions that are used in the simulations and the 
application.

--------------------------------------------- MONTE CARLO ---------------------------------------------
File logit.mlink.R produces the simulation results for the SGMM estimator with missing links at random: 
Figure 1(a) and Tables E.1 and E.2.

File logit.misclassified.R produces the simulation results for the SGMM estimator with misclassified 
links: Figure 1(b) and Tables E.3 and E.4.

File Bayesian.R produces the simulation results for the Bayesian estimator with missing links at random: 
Figure F.1.

File GX_observed_Breza.R produces simulation results for the SGMM estimator with ARD, where GX is 
observed and the network distribution is estimated following Breza et al. (2020): Table H.1.

File GX_unobserved_Breza.R produces simulation results for the SGMM estimator with ARD, where GX is 
unobserved and the network distribution is estimated following Breza et al. (2020): Table H.2.

File GX_observed_Alidaee.R produces simulation results for the SGMM estimator with ARD, where GX is 
observed and the network distribution is estimated following Alidaee et al. (2020): Table H.1.

File GX_unobserved_Alidaee.R produces simulation results for the SGMM estimator with ARD, where GX is 
unobserved and the network distribution is estimated following Alidaee et al. (2020): Table H.2.

--------------------------------------------- ADD HEALTH ---------------------------------------------
File AddHealth.data.R cleans the Inschool AddHealth data and the network data used for the application.

File AddHealth.model.R uses the prepared data from AddHealth.data.R and produces the results for the 
application: Tables G.1, G.2, G.3, and G.4, as well as Figures 2, 3, G.1, G.2, G.3, and G.4.

Information on how to obtain Add Health data files is available on the Add Health website
(www.cpc.unc.edu/addhealth)

Codebook for the Add Health data is available at: https://addhealth.cpc.unc.edu/documentation/codebooks/
