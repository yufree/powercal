# Power Analysis for Metabolic Phenotyping Studies #

Scripts for power analysis 

So far mostly MATLAB code, but python 3.5 versions will be available in the near future...

Main Files and description:

 - Tutorial: Script with examples on how to run the code
 - Pcalc_2Group: Power calculations for difference in mean between 2 groups using ANOVA
 - Pcalc_Continuous: Power calculations for linear regression against a continuous outcome
 - PCalc_Logistic: Power calculations for logistic regression
 - SurfacePlot: Function to plot power analysis results (single variable)
 - TutorialData.mat: Contains a dataset for use in tutorial (C. elegans HRMAS dataset)
 
 Auxiliary files and functions:

 - simulateLogNormal: Function to sample from a multivariate lognormal fitted to the data
 - calcConfMatrixUniv: Confusion matrix calculation function
 - bh_fdr: Function for Benjamini-Hochberg and Benjamini Yekutieli (credits to David Groppe: http://www.mathworks.com/matlabcentral/fileexchange/27418-fdr-bh)

For further information contact:
- gscorreia89_at_gmail.com