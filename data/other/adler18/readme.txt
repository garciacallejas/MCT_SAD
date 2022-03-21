To run the scripts, open CompRegress_wrapper.r and set the 
working directory to your path. You can then run each of the 
component scripts one by one or source it all at once. You
may need to install some packages you don't have, but once 
you do it should run all the way through.

Here is a description of the column headers in the data file,
CompetitionRegressionData071517.csv:

Paper Key #:A unique identifier given to each study
First Author:First author of each study
Year:Calendar year each study published
Data.Extractor:	Name of person who extracted the data for this review
Response variable:How individual plant performance was measured
Greenhouse or field experiment:	Setting of the study
Ecosystem:  Description of vegetation type where the study was conducted
Longterm (>1 season): Study duration	
Manipulated or natural densities?: Source of variation in competitor densities
Focal.Species: Scientific name of the focal species
Focal.Lifestage: Life stage at which focal species' performance measured
Focal.Form: Life form of focal species
Focal.Origin: Whether the focal species is native or non-native
Comp.Species: Scientific name of the competitor species
Comp.Lifestage: Life stage of the competitor species	
Comp.Form: Life form of the competitor
Comp.Origin: Native vs exotic status of the competitor
treatment: Some studies measured competition under a variety of different treatment conditions
nonsig.zero: A value of one indicates that the study assigned coefficient a value of zero if they were not significantly different from zero
competition coefficient: the per capita effect of the competitor on the focal species
Negative coefs mean competition: "Yes" indicates that negative coefficients represent competitive effects, and positive coefficients represent facilitative effects	
standard error of coefficient: standard error of the competition coefficient
Description of Uncertainty Metric: description of metric used to describe uncertainty in the competition coefficient (if any)
P Value: p-value of the competition coefficient (null hypothesis=no different from zero)
DF: degrees of freedom used to estimate competition coefficient
CI lower: Lower bound of confidence interval on the competition coefficient
CI upper: Upper bound of confidence interval on the competition coefficient	
Comments: Any notes about the study
