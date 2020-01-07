# OCP
One-Class Peeling Method


In Phase I of Statistical Process Control (SPC), control charts are often used as outlier detection methods to assess process stability.  Many of these methods require estimation of the covariance matrix, are computationally infeasible, or have not been studied when the dimension of the data, $p$, is large. We propose the one-class peeling (OCP) method, a flexible framework that combines statistical and machine learning methods to detect multiple outliers in multivariate data.  The OCP method can be applied to Phase I of SPC, does not require covariance estimation and is well-suited to high-dimensional datasets with a high percentage of outliers.  Our empirical evaluation suggests that the OCP method performs well in high dimensions and is computationally more efficient and robust than existing methodologies. We motivate and illustrate the use of the OCP method in a Phase I SPC application on a $N=354$, $p = 1,917$ dimensional dataset containing Wikipedia search results for NFL players, teams, coaches, and managers. The example data set and R functions, \textit{OCP.R} and \textit {OCPLimit.R}, to compute the respective OCP distances and thresholds are available in the supplementary materials.
