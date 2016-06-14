# Model for identifying biased SNPs for subsequent removal

## Overview

The code in this directory implements a model for identifying
biased SNPs described in section 1.3.1 of the supplementary 
material. There are two separate implementations of the model,
both of which are run using R and should produce very similar results:
1. The original implementation released with the paper, 
available in the directory [orig](orig). The code in this directory
is updated from that available on the *Genome Research* website to fix
several bugs and improve clarity.
2. An implementation of the model using [JAGS](http://mcmc-jags.sourceforge.net/). JAGS is a general framework
for simulation from Bayesian hierarchical models using MCMC. Here we
call JAGS from R using the rjags package.

## Dependencies

* Original code: the optparse R package
* JAGS implementation: the rjags R package

## Which implementation should you choose? 

The choice is up to you.
The first implementation has been more thoroughly tested but is 
relatively slow and the code is specific to this paper so it may
be difficult to troubleshoot.
The second implementation should be significantly faster and
easier to troubleshoot, but should be considered a preview release 
in that it may contain bugs.
In the longer term development will focus on the JAGS implementation.