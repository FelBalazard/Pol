# Pol
The code contained in this repository is the code used in the article "Policy-aware Evaluation of treatment personalization" by F.Balazard G.Biau P.Ravaud R.Porcher.
It consists of 7 R scripts that have the corresponding subsection at the end of the file name.

pol1_3.1.R is used to show the problems linked with the choice of the plug-in estimate of the optimal policy. It implements our procedure. Figure 1 comes from this script.

scenario1_3.4_4.3.R shows that the max lower bound policy deals with the problems of pol1. The quantities of interest are estimated for several choice of threshold. Figure 3 and 6 come from this script.

scenario2_3.4 is almost the same as the previous script except for the choice of parameter and the focus on the correlation term.

zdelta_3.2.R is used to plot figure 2.

TypeIerror_4.1 is used to plot figure 4 on Type I error control.

coverage_4.2.R is used to plot figure 5 on coverage probability.

multivariate_7.R is used to apply our method and its extension to the multivariate case to the IST trial.
