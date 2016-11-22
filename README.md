# Allele-specific expression

This repository contains code to perform allele-specific expression analysis (ASE) using hierarchical Bayesian model that utilizes the information on both DNA and RNA variation ([Skelly et al. 2011](https://dx.doi.org/10.1101/gr.119784.110)).

This particular instruction is adopted to analyze the data of *Capsella bursa-pastoris* (Brassicaceae). 


## Obtain count data

The allelic expression counts are generated using [ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) with the following steps:


### Generate reference VCF

Extract one sample from a multisample VCF.

```
java -Xmx4g -jar GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R reference.fa \
  -V multisample_SNPs.vcf \
  -sn sample1 \
  -selectType SNP \
  -o sample1.vcf
```

ASEReadCounter counts only REF ALT alleles and can produce incorrect counting for samples that are very divergent from the reference in the case when its both alleles are in ALT field. 
[make_ASEReadCounter_InputVCF.py](make_ASEReadCounter_InputVCF.py) modifies a VCF file by substituting REF ALT fields with two alleles from a heterozygous site. Homozygous sites are skipped.

```
python make_ASEReadCounter_InputVCF.py -i sample1.vcf -o sample1.modif.vcf
```

### Perform allele counting

The allelic expression counts are generated using [ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php):

```
java -Xmx4g -jar GenomeAnalysisTK.jar \
  -T ASEReadCounter \
  -R reference.fa \
  -I sample1.bam \
  -sites sample1.modif.vcf \
  -U ALLOW_N_CIGAR_READS \
  -minDepth 10 \
  --minBaseQuality 20 \
  --minMappingQuality 30 \
  --includeDeletions \
  --outputFormat RTABLE 
  -o sample1.ASEReadCounter.table
```

These allelic counts are then phased using the information on the phasing state:

```
python phaseASEReadCounter_GeneSplit.py -i sample1.ASEReadCounter.table -r reference.file -g genes.bed -o _sample1.output
```

For the format of all input and output files see the code of [phaseASEReadCounter_GeneSplit.py](phaseASEReadCounter_GeneSplit.py). The code to obtain the `reference.file`is provided in the [genome-phasing repository](https://github.com/evodify/genome-phasing).


## Run the scripts by Skelly et al.

The script [phaseASEReadCounter_GeneSplit.py](phaseASEReadCounter_GeneSplit.py) produces several outputs, one of which is used as input for the scripts originally developed by [Skelly et al. (2011)](https://dx.doi.org/10.1101/gr.119784.110). These scripts have been forked from https://github.com/daskelly/ase and slightly modified. The original scripts used for modification are provided in the folder [skelly_original](skelly_original).


### Remove biased SNPs

All the scripts are provided in the folder [find_biased_SNPs](find_biased_SNPs).

#### Remove biased DNA SNPs

The script [find_unbiasedDNA.R](find_biased_SNPs/find_unbiasedDNA.R) takes as input the count data from the genomic DNA where no ASE is expected and removes all SNPs that show highly biased allelic read counts. Diagnostic plots are produced during the output.

To filter out highly biased SNPs run:

```
Rscript find_unbiasedDNA.R DNAinput1.csv DNAinput2.csv DNAinput3.csv 
```

#### Remove biased RNA SNPs

The script [extract_unbiasedRNA.sh](find_biased_SNPs/extract_unbiasedRNA.sh) extracts from the RNA data SNPs that showed no bias in the DNA data.

```
sh extractUnbiased.sh DNAinput1_unbiased.csv RNAinput.csv
``` 

### DNA model

#### Run the DNA model
The DNA counts data is used to estimates overdispersion in the null data where 
should be no ASE.

```
Rscript DNAmodel.R DNAinput1_unbiased.csv DNAinput1_unbiased.gz 100000 50 2000
```
where 100000 is the number of iterations of MCMC, 50 is thin interval, 2000 is the number of scaling iterations. These parameters can be changed.

#### Extract the results of the DNA model

The script [DNAmodel_exctract_results.R](DNA_model/DNAmodel_exctract_results.R) enables to extract the values of a.hat and d.hat that will be needed to run RNA model. It also outputs distributions of a and d to check if there were any problems during the estimate.

```
Rscript DNAmodel_exctract_results.R DNAinput1_unbiased.gz DNAinput2_unbiased.gz DNAinput3_unbiased.gz
```

The values of a.hat and d.hat will be written to the file 'a.hat_d.hat.txt', there also will be a JPEG file with distributions of a.hat and d.hat for each input file. Different runs of DNA model on the same input file should converge to the same values with very similar distributions.

### RNA model

#### Run the RNA model

```
Rscript RNAmodel.R --n.iter=200000 --thin=50 --n.scaling.iter=2000 --a.hat=80 --d.hat=160 --max.rounds.of.scaling=8 RNAinput1_unbiased.csv RNAinput1_unbiased.gz
```
`--a.hat=` and `--d.hat=` use the values obtained in the previous stem for the DNA model.

To see all the possible options, run `Rscript RNAmodel.R  --help` 

#### Extract the results of the RNA model

Extract genes that show posterior probability of ASE higher than 0.99:

```
Rscript RNAmodel_exctract_results.R RNAinput1_unbiased.gz RNAinput2_unbiased.gz RNAinput3_unbiased.gz
```

The threshold 0.99 can be changed by editing `cutoff <- 0.99` in [RNAmodel_exctract_results.R](RNA_model/RNAmodel_exctract_results.R).

Optionally, the results of ASE analysis (files `*_signASE.csv`,`*_allASE.csv`) can be merged with all the expression information (`Genes_*` file generated by [phaseASEReadCounter_GeneSplit.py](phaseASEReadCounter_GeneSplit.py)):

```
sh merge_ASE_Gene.sh RNAinput1_unbiased_allASE.csv Genes_RNAinput1.csv
```

## Possible problems

Both [DNAmodel_exctract_results.R](DNA_model/DNAmodel_exctract_results.R) and [RNAmodel_exctract_results.R](RNA_model/RNAmodel_exctract_results.R) source [readGzippedMcmcOutput.R](DNA_model/readGzippedMcmcOutput.R) that may not work correctly with R/3.0.0 and later. The error message states `In readLines(file.gz): seek on a gzfile connection returned an internal error`. **It worked only in R/2.13.0** for me.
