 # Allele-specific expression

This repository contains code to implement the model of allele-specific 
expression (ASE) described in Skelly et al. 2011 Genome Research 
[doi:10.1101/gr.119784.110](https://dx.doi.org/10.1101/gr.119784.110).

The directory [biased_SNPs](biased_SNPs) provides code to implement the model
for detecting biased SNPs from DNA data. These SNPs are then filtered out for all
subsequent analyses. This model is described in section S1.3.1 of the 
supplementary material. See the README in that directory for more details.

The directory [DNA_model](DNA_model) provides code to implement the model
for genomic DNA read counts that estimates overdispersion in this "null" data where 
no genes should show ASE. This model is described in section 1.3.2 of the 
supplementary material. See the README in that directory for more details.

The directory [RNA_model](RNA_model) provides code to implement the model
to detect ASE in read counts derived from RNA. This model is described in section 1.3.3 
of the supplementary material. See the README in that directory for more details.

[This tutorial](tutorial.pdf) (PDF) provides an overview of how to use
scripts in the `DNA_model/orig` and `RNA_model/orig` directories to implement
our statistical model for ASE. These are the scripts that were originally
released with the paper and are available on the *Genome Research*
website as supplementary information. They have a few corrections from
the exact code published on the *Genome Research* website to account for
bugs discovered after publication as well as changes to dependencies 
that the code utilizes.
