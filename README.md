# Allele-specific expression

This repository contains code to perform allele-specific expression analysis (ASE) using hierarchical Bayesian model that utilizes the information on both DNA and RNA variation ([Skelly et al. 2011](https://dx.doi.org/10.1101/gr.119784.110)).

This particular instruction is adopted to analyze the data of *Capsella bursa-pastoris* (Brassicaceae).

## Map reads to the reference genome

### Mask the variant sites in the reference

Mapping bias is a major problem in this analysis. Masking the variant sites helps to reduce the mapping bias. Here is one of the ways you can approach this problem.

#### Create a BED file of regions to mask

Create a BED file from SNPs and Indels VCFs:
```
zcat SNPs.vcf.gz | grep -v ^# | awk '{printf "%s\t%s\t%s\n", $1,$2-1,$2}' > SNPs.bed

zcat INDELs.vcf.gz | grep -v ^# | awk 'length($4) != 1 {printf "%s\t%s\t%s\t\n", $1,$2-1,$2+length($4)-1}' > heter_INDELs.bed
```
*Be careful with masking indels, if there are many indels, masking them may result is spurious mapping later.*

Using [bedtools](http://bedtools.readthedocs.io/en/latest/), merge the overlapping intervals in the BED files
```
bedtools merge -i SNPs.bed > SNPs.merged.bed
bedtools merge -i INDELs.bed > INDELs.merged.bed
```
Merge the BED files obtained from VCFs and the BED file describing repetitive regions:
```
cat SNPs.merged.bed INDELs.merged.bed repetitiveRegions.bed | sort -V > SNPs-INDELs_repeats.bed
```
Merge the overlapping intervals:
```
bedtools merge -i SNPs-INDELs_repeats.bed > SNPs-INDELs_repeats.merged.bed
```
#### Mask variant sites and repetitive regions in the reference

```
bedtools maskfasta -fi reference.fa -bed SNPs-INDELs_repeats.merged.bed -fo reference_masked.fa
```

### Map the reads

#### Prepare the reference

This analysis will require both normal reference and masked reference files:

```
java -Xmx4g -jar CreateSequenceDictionary.jar R=reference_masked.fa O=reference_masked.dict
samtools faidx reference_masked.fa

java -Xmx4g -jar CreateSequenceDictionary.jar R=reference.fa O=reference.dict
samtools faidx reference.fa
```


#### Map with stampy

[Stampy](http://www.well.ox.ac.uk/project-stampy) requires python2.6.

Prepare the reference files:

```
./stampy.py --species=referenceName --assembly=referenceMasked -G referenceMasked reference_masked.fa.gz
./stampy.py -g referenceMasked -H referenceMasked

```

Perform the alignment of both DNA and RNA reads to the reference.
```
./stampy.py -t16 -g referenceMasked -h referenceMasked -o sample1.sam -M sample1_R1.fastq.gz sample1_R2.fastq.gz
```
*where -t16 means to use 16 cores in a multithreaded processor.*

#### Prepare BAM file

Convert SAM to BAM and sort BAM file with [samtools](http://www.htslib.org/):
```
samtools view -bS -@ 16 sample1.sam > sample1.bam
samtools sort -@ 16 sample1.bam -o sample1_sorted.bam
samtools index sample1_sorted.bam
```

Add groups:
```
java -Xmx6g -jar ~/Programs/picard-tools/AddOrReplaceReadGroups.jar \
  I=sample1_sorted.bam \
  O=sample1_sorted_GroupAdded.bam \
  RGLB=speciesName RGPL=illumina RGPU=DJG63NM1 RGSM=TY118
samtools index sample1_sorted_GroupAdded.bam
```
*RGPU should be specified with lane id (can be seen in a fastq file name or in reads information).*

Mark duplicated reads:
```
java -Xmx6g -jar ~/Programs/picard-tools/MarkDuplicates.jar \
  I=sample1_sorted_GroupAdded.bam \
  O=sample1_sorted_GroupAdded_MarkDupl.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
  M=sample1_sorted_GroupAdded_MarkDupl.metrics
```
Split N cigar reads:
```
java -Xmx6g -jar ~/Programs/GATK/GenomeAnalysisTK.jar \
 -T SplitNCigarReads \
 -R ~/reference/Crubella_183_SNPs-INDELs_repeats_masked.fa  \
 -I sample1_sorted_GroupAdded_MarkDupl.bam \
 -U ALLOW_N_CIGAR_READS \
 -o sample1_sorted_GroupAdded_MarkDupl_split.bam
```
*If you use STAR and some other alignment software, you also need to include mapping quality reassignment*

## Obtain count data

The allelic expression counts are generated using [ASEReadCounter](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) with the following steps:


### Generate reference VCF

Extract one sample from a multisample VCF that is obtained from the DNA sequence data (or RNA data).

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
  -I sample1_sorted_GroupAdded_MarkDupl_split.bam \
  -sites sample1.modif.vcf \
  -U ALLOW_N_CIGAR_READS \
  -minDepth 10 \
  --minBaseQuality 20 \
  --minMappingQuality 30 \
  --includeDeletions \
  --outputFormat RTABLE \
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

To get unbiased counts and allelic ratio per gene, the SNP counts can be merged by genes with [merge_by_Genes.sh](/find_biased_SNPs/merge_by_Genes.sh):

```
sh merge_by_Genes.sh RNAinput1_unbiased.csv
```

To keep only SNPs that overlap in DNA and RNA data, remove DNA SNPs that do not exist on the unbiased RNA data:

```
sh extractUnbiased.sh RNAinput1.csv_unbiased.csv DNAinput1_unbiased.csv
```

The file `DNAinput1_unbiased.csv_unviased.csv` is then used to for the DNA model.


### DNA model

#### Run the DNA model
The DNA counts data is used to estimates overdispersion in the null data where
should be no ASE.

```
Rscript DNAmodel.R DNAinput1_unbiased.csv_unbised.csv DNAinput1_unbiased.gz 200000 100 5000
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
Rscript RNAmodel.R --n.iter=200000 --thin=100 --n.scaling.iter=5000 --a.hat=80 --d.hat=160 --max.rounds.of.scaling=8 RNAinput1_unbiased.csv RNAinput1_unbiased.gz
```
`--a.hat=` and `--d.hat=` use the values obtained in the previous stem for the DNA model.

To see all the possible options, run `Rscript RNAmodel.R  --help`

#### Extract the results of the RNA model

Extract genes that show posterior probability of ASE higher than 0.99:

```
Rscript RNAmodel_exctract_results.R RNAinput1_unbiased.gz RNAinput2_unbiased.gz RNAinput3_unbiased.gz
```

The significance threshold 0.99 can be changed by editing `cutoff <- 0.99` in [RNAmodel_exctract_results.R](RNA_model/RNAmodel_exctract_results.R). The corresponding FDR threshold is reported in the last line of the output file.
This script also report the distribution of estimates statistics which should be explored for the convergence between different runs.


Optionally, the script [merge_ASE_with_Gene.sh](RNA_model/merge_ASE_with_Gene.sh) can be used to merge the results of ASE analysis (files `*_signASE.csv`,`*_allASE.csv`) with all the expression information (`Genes_*_unviased.csv` files):

```
sh merge_ASE_with_Gene.sh RNAinput1_unbiased_allASE.csv Genes_RNAinput1_unbiased.csv
```

## Possible problems

Both [DNAmodel_exctract_results.R](DNA_model/DNAmodel_exctract_results.R) and [RNAmodel_exctract_results.R](RNA_model/RNAmodel_exctract_results.R) source [readGzippedMcmcOutput.R](DNA_model/readGzippedMcmcOutput.R) that may not work correctly with R/3.0.0 and later. The error message states `In readLines(file.gz): seek on a gzfile connection returned an internal error`. It worked only in **R/2.13.0** for me.
Another workaround is to use uncompressed input (gunzip), but [readGzippedMcmcOutput.R](DNA_model/readGzippedMcmcOutput.R) has to be edited for this type of input.

**DISCLAIMER:** USE THESE SCRIPTS AT YOUR OWN RISK. I MAKE NO WARRANTIES THAT THESE SCRIPTS ARE BUG-FREE, COMPLETE, AND UP-TO-DATE. I AM NOT LIABLE FOR ANY LOSSES IN CONNECTION WITH THE USE OF THESE SCRIPTS.
