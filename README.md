# Allele-specific expression

This repository contains code to perform allele-specific expression analysis (ASE) using hierarchical Bayesian model that utilize the information on both DNA and RNA variation and allow both global and locus-specific inferences about allele-specific expression. This particular instruction is adopted to analyze the data of *Capsella bursa-pastoris* (Brassicaceae). 

## Obtain count data

The allelic expression counts are generated using [ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php):

```
java -Xmx4g -jar GenomeAnalysisTK.jar \
  -T ASEReadCounter \
  -R reference.fa \
  -I input.bam \
  -sites input_SNPs.vcf \
  -U ALLOW_N_CIGAR_READS \
  -minDepth 10 \
  --minBaseQuality 20 \
  --minMappingQuality 30 \
  --includeDeletions \
  --outputFormat RTABLE 
  -o ASEReadCounter.table
```

These allelic accounts are then phased using the information on the phasing state:

```python phaseASEReadCounter_GeneSplit.py -i ASEReadCounter.table -r reference.file -g genes.bed -o _output
```

The code to obtain the `reference.file`is provided in the [genome-phasing repository](https://github.com/evodify/genome-phasing). For the format of all input and output files see the code of [phaseASEReadCounter_GeneSplit.py](phaseASEReadCounter_GeneSplit.py).

The script *phaseASEReadCounter_GeneSplit.py* produces several outputs, one of which is used as input to the scripts originally developed by Skelly et al. 2011. Genome Research [doi:10.1101/gr.119784.110](https://dx.doi.org/10.1101/gr.119784.110). These scripts have been forked from [daskelly/ase](https://github.com/daskelly/ase) and slightly modified.

## Run the scripts by Skelly et al.


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
