% Tutorial for power analysis codes 
%
% This tutorial has two parts. 
% Part I - Test the code for a single combination of parameters (sample size, effect size etc), to check the code runs correctly.
% 	This typically runs in less than a minute on a typical 4-core desktop.
% Part II - Run the code for a grid of many sample/effect sizes and generate a surface plot of false negative rates.
%
% 2) single variable at many grid points, aiming to take 1 minute ish... clearly state that to do all variables in the data takes longer!
% Give brief advice on how to speed up e.g. reduce number of variables / repeats / grid points...

%load C. elegans dataset 
load('TutorialData.mat');

% Part I
% Run code for a single effect and sample size combination as a test
effectSizes = [0.5];
sampleSizes = [200];
% give numberreps a number
numberreps = 10;

[diffgroups] = PCalc_2Group(XSRV,effectSizes, sampleSizes, 0.05, 5000, numberreps);
[linearregression] = PCalc_Continuous(XSRV,effectSizes, sampleSizes, 0.05, 5000, numberreps);

% Part II 
% Define a grid of effect sizes and sample sizes to test
effectSizes = [0.05, 0.1, 0.15,0.2, 0.25, 0.3, 0.35];
sampleSizes = [50, 100, 200, 250, 350, 500, 750, 1000];
numberreps= 10;

% Calculat for a subset of 4 variables (less than 20 seconds on 4-core desktop for each analysis)
[diffgroups] = PCalc_2Group(XSRV(:, 1:4),effectSizes, sampleSizes, 0.05, 5000, numberreps);
[linearregression] = PCalc_Continuous(XSRV(:, 1:4),effectSizes, sampleSizes, 0.05, 5000, numberreps);

% Surface plot function (see details in bottom of tutorial)
SurfacePlot(diffgroups, 2, 4,2 , sampleSizes, effectSizes,numberreps)


% Run the code for all variables. Each analysis takes around 1h on a 4 core desktop. To speed up, use less effect and sample 
% sample sizes and a smaller number of repeats
[diffgroups] = PCalc_2Group(XSRV,effectSizes, sampleSizes, 0.05, 5000, numberreps);
[linearregression] = PCalc_Continuous(XSRV,effectSizes, sampleSizes, 0.05, 5000, numberreps);

%{
Using the SurfacePlot function to visualize results 
SurfacePlot(output, variable,metric,correction, sizeeff,samplsizes,nreps)
Output is the structure returned from the simulator, variable is the index of variable to plot
metric is the to display and correction the type of multiple testing correction to 
visualize.

Metric options:
1 - True positive Rate
2 - False Positive Rate
3 - True Negative Rate
4 - False Negative Rate
Correction:
1 - No correction
2 - Bonferroni
3 - Benjamini-Hochberg
4 - Benjamini-Yekutieli

The example line below will open the False Negative Rate surface for
variable number 2 without multiple testing correction
%}
SurfacePlot(diffgroups, 2, 4,2 , sampleSizes, effectSizes,numberreps)

