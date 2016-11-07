# Model for identifying biased SNPs for subsequent removal

## Overview

The code in this directory implements a model for identifying
biased SNPs described in section 1.3.1 of the supplementary 
material in 
[doi:10.1101/gr.119784.110](https://dx.doi.org/10.1101/gr.119784.110). 
There are two separate implementations of the model,
both of which are run using R and should produce very similar results:

1. The original implementation used in the paper, 
available in the directory [orig](orig). The code in this directory
was used to filter out putatively biased SNPs in that analysis.
2. An implementation of the model using [JAGS](http://mcmc-jags.sourceforge.net/). JAGS is a general framework
for simulation from Bayesian hierarchical models using MCMC. Here we
call JAGS from R using the `runjags` and `rjags` packages.

## Dependencies

* [orig](orig): the `optparse` R package
* [JAGS](JAGS): the `rjags` and `runjags` R packages

## Which implementation should you choose? 

The choice is up to you.
The first implementation has been more thoroughly tested but is 
relatively slow and the code is specific to this paper so it may
be difficult to troubleshoot.
The second implementation should be significantly faster and much
easier to troubleshoot, but should be considered a preview release 
in that it may still contain bugs.

Over the longer term, development will focus on the JAGS implementation
and the original code will be archived.

## Note about the model

This model is not strictly necessary for running the model of ASE
described in the paper. Nevertheless, it is important to make sure that 
SNPs that show highly biased allelic read counts in genomic DNA are
removed. In our study, we expected that allelic read counts would show 
some variability beyond that expected from statistical sampling, 
due to factors such as the many steps involved in preparation of sequencing
libraries. Read counts at a minority of SNPs in our genomic DNA data 
appeared highly biased, even after removing SNPs where reads generated 
*in silico* did not show 50/50 mapping of alleles.
The model above was thus used to filter out such highly biased SNPs.

The results of this model should be fairly similar to the simpler 
strategy of filtering out SNPs that have small *p*-values according
to the binomial exact test. In our experience, 
the model implemented here is slightly
more conservative (i.e. more SNPs are called as biased) than 
removing those SNPs that have a Bonferroni-corrected *p*-value
less than 0.05 according to the binomial exact test.

## Example code

Here is some code to run a test analysis from within R. 
Using the original model:
```
setwd("orig")
source("run_test.R")
source("../calculate_posterior_prob_biased.R")
```

Using the JAGS model:
```
setwd("rjags")
source("run_test.R")
source("../calculate_posterior_prob_biased.R")
```

Note that these simulations use "toy" data simulated under the
same model as that being inferred, so they do not provide
any indication of whether *your* data fits the model.

For a real analysis of non-simulated data,
you could start from either of the `run_test.R` scripts above
and modify it to read in your data and run the model.
This should be straightforward as the `run_test.R` scripts are less than 
25 lines of code. Then, as above, you can use code in
`calculate_posterior_prob_biased.R` to calculate the posterior probability
that each SNP is biased. You will need to change 
the assignment of `Y` and `N` in lines 22-23 to
reflect your data rather than the simulated data.
Then you will filter out SNPs with high posterior probability of being
biased from further consideration.
