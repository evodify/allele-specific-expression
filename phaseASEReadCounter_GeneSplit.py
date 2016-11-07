#!/usr/bin/env python2

"""
This script phases allele specific expression data produced by ASEReadCounter (GATK).
A reference file of phased homeologues and gene intervals in BED format are required.

# command:
python phaseASEReadCounter_GeneSplit.py -i ASEReadCounter.table -r reference.file -g genes.bed -o _output

# contact:
Dmytro Kryvokhyzha <dmytro.kryvokhyzha@evobio.eu>


#ASEReadCounter.table:
contig  position    variantID   refAllele   altAllele   refCount    altCount    totalCount  lowMAPQDepth    lowBaseQDepth   rawDepth    otherBases  improperPairs
scaffold_1  1904    .   C   A   26  11  37  0   0   39  0   2
scaffold_1  1919    .   C   A   18  17  35  0   0   36  0   1
scaffold_1  2556    .   T   G   11  12  23  0   1   26  0   2
scaffold_1  2594    .   A   T   17  2   19  0   1   20  0   0
scaffold_1  2618    .   A   T   4   14  18  0   0   18  0   0
scaffold_1  2629    .   T   C   5   14  19  0   0   19  0   0
scaffold_1  2639    .   A   C   13  5   18  0   0   18  0   0
scaffold_1  2686    .   T   C   3   13  16  0   0   16  0   0
scaffold_1  2825    .   T   A   18  22  40  0   0   43  0   3

# reference.file
CHROM  POS sample1.SNPs.haplotype_A   sample1.SNPs.haplotype_B
scaffold_1  1898    T   T
scaffold_1  1919    A   C
scaffold_1  2556    G   T
scaffold_1  2594    N   N
scaffold_1  2618    T   A
scaffold_1  2629    C   T
scaffold_1  2639    A   C
scaffold_1  2686    C   T
scaffold_1  2825    A   T

# genes.bed
scaffold_1  1910    2000    Carubv10010485m.g
scaffold_1  2600    3000    Carubv10011809m.g
scaffold_1  7595    9248    Carubv10011354m.g

# SNPs_output
contig  position    variantID   HomeologueA HomeologueB HomeologueACount    HomeologueBCount    TEST_DNA_stampy025.csv.Ratio
scaffold_1  1919    .   A   C   17  18  0.485714285714
scaffold_1  2618    .   T   A   14  4   0.777777777778
scaffold_1  2629    .   C   T   14  5   0.736842105263
scaffold_1  2639    .   A   C   13  5   0.722222222222
scaffold_1  2686    .   C   T   13  3   0.8125
scaffold_1  2825    .   A   T   22  18  0.55

# Genes_output
Gene    numberSNPS  HomeologueAmeanCount    HomeologueBmeanCount    TEST_DNA_stampy025.csv.meanRatio
Carubv10010485m.g   1   17.0    18.0    0.485714285714
Carubv10011809m.g   5   15.2    7.0 0.719868421053

# SNPs_OutOfGenes_output
contig  position    variantID   HomeologueA HomeologueB HomeologueACount    HomeologueBCount    TEST_DNA_stampy025.csv.Ratio
scaffold_1  1904    .   A   C   11  26  0.297297297297
scaffold_1  2556    .   G   T   12  11  0.521739130435

"""

############################ modules ##############################

import argparse, sys
import matplotlib.pyplot as plt
import numpy as np

############################ options ##############################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input_to_phase', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-r', '--phasing_reference', help = 'file containing list of alleles for parental species', type=str, required=True)
parser.add_argument('-g', '--gene_intervals', help = 'file containing gene intervals', type=str, required=True)
args = parser.parse_args()

############################ functions ###########################

def phase_state(GT, RefGT):
  """ define a phasing state relative to the reference """
  if (GT[0] == RefGT[0]) and (GT[1] == RefGT[1]):
    return 'same'
  elif (GT[0] == RefGT[1]) and (GT[1] == RefGT[0]):
    return 'reverse'
  else:
    return 'NA'

############################ script ##############################

RefFile = open(args.phasing_reference, "r")
Refheader = RefFile.readline()
Refwords = RefFile.readline().split()
refChr = Refwords[0].split('_')[1]
refPos = Refwords[1]


geneIntervals = open(args.gene_intervals, "r")
geneIntwords = geneIntervals.readline().split()
geneIntChr = geneIntwords[0].split('_')[1]
geneIntStart = geneIntwords[1]
geneIntEnd = geneIntwords[2]
geneName =  geneIntwords[3]

counter = 0
counterGenes = 0
counterGenesOut = 0
counterGenesWithin = 0

homeologue1Count = []
homeologue2Count = []
ratioList = []
ratioAll = []

# create output files
outputSNPs = open("SNPs" + args.output, 'w') # phased SNPs with counts and allelic ratio
outputGenes = open("Genes" + args.output, 'w') # average counts and allelic ratio per gene
skelly =  open("Skelly" + args.output, 'w') # input for the scripts by Skelly et al. 2011
outputSNPsOut = open("SNPs_OutOfGenes" + args.output, 'w') # phased SNPs that are outside of genes' coordinates

outputSNPs.write("contig\tposition\tvariantID\tHomeologueA\tHomeologueB\tHomeologueACount\tHomeologueBCount\t%s.Ratio\n" % args.input_to_phase)
skelly.write("gene\tvariantID\tHomeologueAcount\tHomeologueBcount\n")
outputSNPsOut.write("contig\tposition\tvariantID\tHomeologueA\tHomeologueB\tHomeologueACount\tHomeologueBCount\t%s.Ratio\n" % args.input_to_phase)
outputGenes.write("Gene\tnumberSNPS\tHomeologueAmeanCount\tHomeologueBmeanCount\t%s.meanRatio\n" % args.input_to_phase)

with open(args.input_to_phase) as datafile:

  # make header line:
  line1 = datafile.readline()
  print('Phasing ...')
  
  # read the file line by line:
  for line in datafile:
    words = line.split()
    Chr = words[0].split('_')[1]
    Pos = words[1]
    chr_pos_id = words[0:3]
    GT = words[3:5]
    GTTotalCount = words[5:8]
    rowName = '\t'.join(str(e) for e in chr_pos_id)

    # find corresponding gene coordinates if such exist
    while int(geneIntChr) < int(Chr) or (int(geneIntChr) == int(Chr) and int(geneIntEnd) < int(Pos)):
      geneIntwords = geneIntervals.readline().split()
      geneIntChr = geneIntwords[0].split('_')[1]
      geneIntStart = geneIntwords[1]
      geneIntEnd = geneIntwords[2]

    # find corresponding reference genotypes if such exist
    while int(refChr) < int(Chr) or (int(refChr) == int(Chr) and int(refPos) < int(Pos)):
      Refwords = RefFile.readline().split()
      refChr = Refwords[0].split('_')[1]
      refPos = Refwords[1]
    if int(refChr) == int(Chr) and int(refPos) == int(Pos):
      RefGT = Refwords[2:4]
    else:
      RefGT = ['N', 'N']
      
    # define the phasing state
    phasedS = phase_state(GT, RefGT)
    
    #print Chr, Pos, GT, RefGT, phasedS # for debugging
    
    if phasedS == 'same':
      homeologueGT = '\t'.join(str(e) for e in GT)
      GTTotalCountR = GTTotalCount
      homeologueCount = '\t'.join(str(e) for e in GTTotalCount[0:2])
      ratio = float(GTTotalCount[0])/(float(GTTotalCount[0])+float(GTTotalCount[1]))
    elif phasedS == 'reverse':
      GTR = [GT[1], GT[0]]
      GTTotalCountR = [GTTotalCount[1], GTTotalCount[0]]
      homeologueGT = '\t'.join(str(e) for e in GTR)
      homeologueCount = '\t'.join(str(e) for e in GTTotalCountR)
      ratio = float(GTTotalCount[1])/(float(GTTotalCount[0])+float(GTTotalCount[1]))
    else:
      GTR = ['N', 'N']
      GTTotalCountR = ['NA','NA']
      homeologueGT = '\t'.join(str(e) for e in GTR)
      homeologueCount = '\t'.join(str(e) for e in GTTotalCountR)
      ratio = 'NA'
      continue # do not process missing data
    
    # process SNPs
    if int(geneIntChr) == int(Chr) and int(geneIntStart) < int(Pos) and int(geneIntEnd) > int(Pos):
      # if within a gene
      geneNameOverlap = geneIntwords[3]
      counterGenesWithin += 1
      outputSNPs.write("%s\t%s\t%s\t%s\n" % (rowName, homeologueGT, homeologueCount, ratio))
      row = [geneNameOverlap, Chr+'_'+Pos, homeologueCount]
      rowP = '\t'.join(str(e) for e in row)
      skelly.write("%s\n" % rowP)
      ratioAll.append(ratio)
    else:
      # if outside of a gene
      geneNameOverlap = "NA"
      counterGenesOut += 1
      outputSNPsOut.write("%s\t%s\t%s\t%s\n" % (rowName, homeologueGT, homeologueCount, ratio)) 

    #print geneIntChr, geneIntStart, geneIntEnd, geneNameOverlap  # for debugging

    # process genes
    if geneNameOverlap == geneName:
      # check if a SNP belongs to the same gene as a previous SNP 
      if ratio != "NA":
        homeologue1Count.append(float(GTTotalCountR[0]))
        homeologue2Count.append(float(GTTotalCountR[1]))
        ratioList.append(float(ratio))
    else:
      # if a SNP belongs to a different gene
      if geneNameOverlap == 'NA':
        # skip intergenic SNPs
        continue
      else:
        if ratioList:
          homeologue1mean = sum(homeologue1Count) / float(len(homeologue1Count))
          homeologue2mean = sum(homeologue2Count) / float(len(homeologue2Count))
          ratioMean = sum(ratioList) / float(len(ratioList))
        else:
          homeologue1mean = homeologueGT[0]
          homeologue2mean = homeologueCount[1]
          ratioMean = ratio
      
      counterGenes += 1
      outputGenes.write("%s\t%s\t%s\t%s\t%s\n" % (geneName, len(homeologue2Count), homeologue1mean, homeologue2mean, ratioMean))
      geneName = geneNameOverlap
      if ratio != "NA":
        homeologue1Count = [float(GTTotalCountR[0])]
        homeologue2Count = [float(GTTotalCountR[1])]
        ratioList = [float(ratio)]

    # track the progress
    counter += 1
    if counter % 100000 == 0:
      print str(counter), "lines processed"
  
  # write the last gene info
  homeologue1mean = sum(homeologue1Count) / float(len(homeologue1Count))
  homeologue2mean = sum(homeologue2Count) / float(len(homeologue2Count))
  ratioMean = sum(ratioList) / float(len(ratioList))
  counterGenes += 1
  outputGenes.write("%s\t%s\t%s\t%s\t%s\n" % (geneName, len(homeologue2Count), homeologue1mean, homeologue2mean, ratioMean))
  geneName = geneNameOverlap
  
  # print some stats
  print "done!\n"
  print str(counterGenesOut), "SNPs outside of genes"
  print str(counterGenesWithin), "SNPs within genes"
  print str(counterGenes), "genes\n"

# Plot a histogram
plt.hist(ratioAll, color="grey", bins=np.arange(0,1.04,0.04))
plt.xlim(0,1)
plt.xticks(np.arange(0,1.1,0.1))
plt.title("Allelic ratio", size = 20)
plt.xlabel("Proportion of the homeologue A")
plt.ylabel("Frequency")
plt.savefig("SNPs" + args.output +".pdf", dpi=90)

outputSNPs.close()
skelly.close()
outputSNPsOut.close()
outputGenes.close()
datafile.close()
geneIntervals.close()
RefFile.close()
