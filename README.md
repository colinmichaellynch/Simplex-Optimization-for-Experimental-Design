# Simplex-Optimization-for-Experimental-Design
I use Sequential Simplex Optimization to Estimate Sample Size for Experiments with a Random Effect 

## Table of Contents 

* [Data](https://github.com/colinmichaellynch/Simplex-Optimization-for-Experimental-Design/blob/main/TrueResponseSurfaceTest.mat)

* [Script](https://github.com/colinmichaellynch/Simplex-Optimization-for-Experimental-Design/blob/main/simplexSearch.m)

* [Video of Simplex in Action](https://github.com/colinmichaellynch/Simplex-Optimization-for-Experimental-Design/blob/main/simplexSearchDesirability.mp4)

## Background 

In this [repository](https://github.com/colinmichaellynch/Sampling-Across-vs-Within-Random-Effects), I discuss a method I developed to find optimal experimental designs when there is a random effect. 'Optimal' in this context means that a design will minimize effort while maximizing balanced accuracy. The main problem with productizing this method is that it is computationally expensive and it does not scale well. It requires running simulations for all combinations of sampling strategies up to some arbitrary limit. As these combinations are 2-dimensional (number of samples within a random effect vs number of samples across a random effect), the total number simulations will be the square of the limit, and so computation time increases exponentially with a larger search area. As I intend on designing an R package centered on implementing this package, this is unacceptable, as many of my users may not have access to powerful computers. My goal here is to find an search method which will find an optimal strategy without so many trials. The first optimization strategy I test is the sequential simplex method. 

## Methods

* I first create a reference response surface where I calculate balanced accuracy for all combinations of sampling strategies. That is, I calculate balanced accuracy for the number of samples within a random effect (the number of within-colony replicates) and across random effects (the number of colonies samples). I also find the effort level for all of these strategies. 

* I set a maximum effort level of 400. All strategies above this threshold are not considered. 

* I find the global optimum. This is the strategy that has the highest balanced accuracy that has an effort level below 100. 

* I transform the response surface with a desirability function. Every strategy above the 400 threshold has a desirablity level set to 0. The global maximum is set to 1. All other strategies have a value between 0 and 1. 

* I use sequential simplex optimization to find the global maximum. This means that I create a random triangle (or simplex in higher dimensions) in the allowed space of strategies. I then run strategies at the corners of this triangle to calculate desireability for these strategies. I then use a set of rules to transform the triangle so that it moves up the desirability gradient. 

* In this video, the x-axis gives the number of ants samples in a colony, the y-axis gives the number of colonies, color shows desirability of a sampling strategy, the bold dotted line shows the maximum effort level, and the black dot is the globabl maximum that the alorithm is searching for. 

https://user-images.githubusercontent.com/61156429/211679683-71e7a3a1-215e-495a-82bd-27eadf509f52.mp4

## Results

* After running this algorithm with 1,000 random starting triangles, the final simplex method gives a balanced accuracy that is within 5% of the global maximum 98.4% of the time with only 1.1% of the number of simulations.  

* This is a much more managable computational load than the exhaustive search, and can be implemented in an R package for selecting optimal sampling strategies.

<p align="center">
  <img src=/Images/simplexPerformance.png>
</p>

## Acknowledgements

I would like to thank my collaborators Dr. Douglas Montgomery and Dr. Daniel McCarville for their valuable input on this project. 
