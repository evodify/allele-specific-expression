#!/usr/bin/env python2

'''
This script creates an input VCF for ASEReadCounter (GATK). ASEReadCounter counts only REF ALT alleles and can produce incorrect counting for samples that are very divergent from the reference and its both alleles are in ALT field. 
This script modifies a VCF file by substituting REF ALT fields with two alleles from a heterozygous site. Homozygous sites are skipped.

# input:

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  sample1
scaffold_1  70699   .   G   A,C 49840.25    PASS    AC=1,1;AF=0.500,0.500;AN=2;BaseQRankSum=-0.190;ClippingRankSum=0.715;DP=78;ExcessHet=0.4885;FS=0.000;InbreedingCoeff=0.3302;MQ=91.67;MQRankSum=-4.218;QD=28.93;ReadPosRankSum=-0.284;SOR=0.738  GT:AD:DP:GQ:PL  1/2:1,39,38:78:99:2678,1582,1469,1186,0,1622
scaffold_1  70866   .   C   T   5891.13 PASS    AC=0;AF=0.00;AN=2;BaseQRankSum=-0.413;ClippingRankSum=-0.518;DP=43;ExcessHet=9.4645;FS=0.000;InbreedingCoeff=-0.2973;MQ=95.67;MQRankSum=1.456;QD=14.88;ReadPosRankSum=0.478;SOR=0.731   GT:AD:DP:GQ:PL  0/0:43,0:43:99:0,130,1622
scaffold_1  70945   .   C   T   13909.41    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=0.970;ClippingRankSum=-0.842;DP=46;ExcessHet=46.7498;FS=0.528;InbreedingCoeff=-0.8462;MQ=92.53;MQRankSum=3.470;QD=14.31;ReadPosRankSum=3.186;SOR=0.617  GT:AD:DP:GQ:PL  0/1:25,21:46:99:735,0,675
scaffold_1  71062   .   G   A   13021.41    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-0.172;ClippingRankSum=0.099;DP=50;ExcessHet=46.7498;FS=1.767;InbreedingCoeff=-0.8462;MQ=92.73;MQRankSum=1.705;QD=15.41;ReadPosRankSum=-0.515;SOR=0.550 GT:AD:DP:GQ:PL  0/1:20,30:50:99:1030,0,579
scaffold_1  71459   .   G   A   45955.25    PASS    AC=2;AF=1.00;AN=2;BaseQRankSum=-2.894;ClippingRankSum=-0.301;DP=47;ExcessHet=0.4885;FS=4.540;InbreedingCoeff=0.3302;MQ=97.68;MQRankSum=0.225;QD=27.70;ReadPosRankSum=9.187;SOR=0.407    GT:AD:DP:GQ:PL  1/1:1,46:47:99:2033,111,0
scaffold_1  71480   .   A   G   49603.25    PASS    AC=2;AF=1.00;AN=2;BaseQRankSum=5.098;ClippingRankSum=-0.761;DP=57;ExcessHet=0.4885;FS=1.641;InbreedingCoeff=0.3302;MQ=97.44;MQRankSum=1.147;QD=35.76;ReadPosRankSum=3.273;SOR=0.489 GT:AD:DP:GQ:PL  1/1:0,57:57:99:2347,171,0
scaffold_1  71548   .   G   A   19906.67    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=1.095;ClippingRankSum=-0.648;DP=65;ExcessHet=29.5512;FS=0.000;InbreedingCoeff=-0.6667;MQ=96.82;MQRankSum=-0.273;QD=15.96;ReadPosRankSum=0.976;SOR=0.683 GT:AD:DP:GQ:PL  0/1:34,31:65:99:977,0,1118

#output:

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  sample1
scaffold_1  70699   .   A   C   49840.25    PASS    AC=1,1;AF=0.500,0.500;AN=2;BaseQRankSum=-0.190;ClippingRankSum=0.715;DP=78;ExcessHet=0.4885;FS=0.000;InbreedingCoeff=0.3302;MQ=91.67;MQRankSum=-4.218;QD=28.93;ReadPosRankSum=-0.284;SOR=0.738  GT:AD:DP    0/1:39,38:78
scaffold_1  70945   .   C   T   13909.41    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=0.970;ClippingRankSum=-0.842;DP=46;ExcessHet=46.7498;FS=0.528;InbreedingCoeff=-0.8462;MQ=92.53;MQRankSum=3.470;QD=14.31;ReadPosRankSum=3.186;SOR=0.617  GT:AD:DP    0/1:25,21:46
scaffold_1  71062   .   G   A   13021.41    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-0.172;ClippingRankSum=0.099;DP=50;ExcessHet=46.7498;FS=1.767;InbreedingCoeff=-0.8462;MQ=92.73;MQRankSum=1.705;QD=15.41;ReadPosRankSum=-0.515;SOR=0.550 GT:AD:DP    0/1:20,30:50
scaffold_1  71548   .   G   A   19906.67    PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=1.095;ClippingRankSum=-0.648;DP=65;ExcessHet=29.5512;FS=0.000;InbreedingCoeff=-0.6667;MQ=96.82;MQRankSum=-0.273;QD=15.96;ReadPosRankSum=0.976;SOR=0.683 GT:AD:DP    0/1:34,31:65

# command:

$ python make_ASEReadCounter_InputVCF.py -i sample1.vcf -o sample1.modif.vcf

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################ modules ##############################

import argparse, sys

############################# options #############################

class CommandLineParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program  #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  print('Creating the output file...')
  fileoutput = open(args.output, 'w')
  
  for line in datafile:
    if line.startswith("#"):  # define header
      fileoutput.write(line)
    else:
      words = line.split()
      GTfield = words[9]
      GTwords = GTfield.split(':')
      GT = GTwords[0].split('/')

      if GT[0] == GT[1]:  # skip if homozygous
        continue
      else:  # extract genotypes for heterozygous sites
        chr_pos = words[0:3]
        ref = [words[3]]
        alt = words[4].split(',')
        refAlt = ref + alt
        inter = words[5:9]
        AD = GTwords[1].split(',')
        DP = GTwords[2]
        ref = refAlt[int(GT[0])]
        alt = refAlt[int(GT[1])]
 
        # process annotation
        if 'AD' not in inter[3]:  # if AD and DP are not present in annotation
          AD = 0
          DP = 0
        else:
          AD = ','.join(str(e) for e in [AD[int(GT[0])]] + [AD[int(GT[1])]])

        # output the results
        chr_posP = '\t'.join(str(e) for e in chr_pos)
        refAltP = '\t'.join(str(e) for e in ref+alt)
        interP = '\t'.join(str(e) for e in inter[0:3])
        fileoutput.write('%s\t%s\t%s\tGT:AD:DP\t0/1:%s:%s\n'  % (chr_posP, refAltP, interP, AD, DP))

    # track the progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
