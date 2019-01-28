 #!/usr/bin/python2

"""
This script generates a custom hybrid reference from a normal reference file and a genotype file of hybrid species. Such reference file can be used for stringent competitive mapping to two genomes simultaneously and subsequent homeologue-specific gene expression analysis.

# input.tab:

#CHR    POS REF sample1 sample2 sample3 sample4 sample5
scaffold_2  1   G   T   N   N   N   N
scaffold_2  2   A   A   A   A   G   A
scaffold_2  3   C   T   C   C   C   C
scaffold_2  5   A   T   T   T   T   T
scaffold_4  1   T   T   T   C   C   T
scaffold_4  2   T   N   T   T   T   T
scaffold_4  3   A   N   A   A   A   A
scaffold_4  4   G   N   G   G   G   G
scaffold_4  6   A   A   A   A   A   A
scaffold_4  7   G   G   G   G   G   G
scaffold_4  8   T   T   T   T   T   T
scaffold_4  9   A   N   A   A   A   A
scaffold_4  10  A   N   A   C   C   C


# input.fasta:

>scaffold_1
TTA
>scaffold_2
GACTA
>scaffold_3
TAGGACTATATTAACCTTAGT
>scaffold_4
TTAGAAGTAA
>scaffold_5
TAGGACTATAACC
>scaffold_6
TAGGACTATAACCTTAGTAGTAACCTGCCGACTATTAGTAGTAACCTGCC
>scaffold_7
TAGGACTATAACCGACTATTAGTAGTAACCTGCC


# output.fasta:

>scaffold_1
TTA
>scaffold_2.1
TATTT
>scaffold_2.2
GACTT
>scaffold_3
TAGGACTATATTAACCTTAGT
>scaffold_4.1
TTAGAAGTAA
>scaffold_4.2
CTAGAAGTAC
>scaffold_5
TAGGACTATAACC
>scaffold_6
TAGGACTATAACCTTAGTAGTAACCTGCCGACTATTAGTAGTAACCTGCC
>scaffold_7
TAGGACTATAACCGACTATTAGTAGTAACCTGCC


# command:

$ python2 hybrid_reference_maker.py -c test.tab -f test-ref.fasta -o output.fasta -s "group1[sample1,sample2];group2[sample3,sample4,sample5]"

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module
import re, random, sys

############################# classes and functions #############################

class FastaParser(object):
    def __init__(self, filename):
        self.filename = filename
        self.name = []
        self.num = []
        self.sequence = []
        self.count = 0

        # read fasta file
        fasFile = open(self.filename, 'r')
        seq = ''
        for line in fasFile:
            words = line.split()
            if words[0].startswith('>'):  # find sequence names
                self.count += 1
                seqName = words[0][1:]
                seqNum = words[0][1:].split("_")[1]
                self.num.append(seqNum)
                self.name.append(seqName)  # append sequence names
                if seq:
                    self.sequence.append(seq)
                seq = ''
            else:
                seq += words[0]  # append sequences
        self.sequence.append(seq)  # append the last sequence.

    #def __len__(self):
        #"""Enables length property for sequences in a file """
        #return self.count

    def __getitem__(self, i):
        """ Enables iteration through the chromosomes/scaffolds by index and names"""
        if isinstance(i, int):  # if index
            if i > self.count or i < 0:
                raise IndexError("Index is out of range")
            return self.sequence[i]
        else:   # if name
            if i not in self.name:
                raise KeyError("No sequence with name %s", i)
            seqIndex = self.name.index(i)
            return self.sequence[seqIndex]


############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-f', '--fastaReference', help = 'name of fasta reference file', type=str, required=True)
parser.add_argument('-c', '--callsFile', help = 'name of genotype calls file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to make a reference for. Specify in the format "group1[sample1,sample2];group2[sample3,sample4,sample5]"', type=str, required=True)
args = parser.parse_args()

# check and append group names and samples
groupNames = args.samples
groups = []
for groupi in groupNames.strip("\"").split(";"):
  groupNames = groupi.split("[")[0]
  groupSample = re.split("\[|\]", groupi)[1]
  groups.append(groupNames)
  vars()[groupNames + "samples"] = calls.checkSampleNames(groupSample,args.callsFile)

############################ script ##############################

fileoutput = open(args.output, 'w')
fastaRef = FastaParser(args.fastaReference)

# to track lines in the input files
fastaNum = 0
CHRprev = 0
CHRnum = 0
fastaNumAdd = 0

print('Opening the file...')

with open(args.callsFile) as callsFile:
  header_words = callsFile.readline().split()

  # index samples according to the header
  for groupName in groups:
    vars()[groupName + "Index"]  = calls.indexSamples(vars()[groupName + "samples"], header_words)

  for line in callsFile:
    words = line.split()
    CHR = words[0]
    POS = int(words[1])
    # prepare the genotypes
    for popName in groups:
      sGT = calls.selectSamples(vars()[popName + "Index"], words) # select genotypes per group
      sGTnoN = [i for i in sGT if i != 'N'] # remove missing data
      if sGTnoN != []:  # if not all GT are missing
        sGTnoNset = list(set(sGTnoN)) # find set of alleles
        random.shuffle(sGTnoNset) # shuffle list to deal with a tie case.
        mostFreqGT = max(sGTnoNset, key=sGTnoN.count) # find the most frequent allele. A tie is solved randomly.
      else:
        mostFreqGT = 'NA'

      if CHR != CHRprev:
        try:
          fastaSeqP = ''.join(str(w) for w in vars()[popName + "fastaSeq"])
          fastaNameP = fastaRef.name[CHRnum-1]
          fileoutput.write(">%s.%s\n%s\n" % (fastaNameP, popName[-1], fastaSeqP))
          fastaNumAdd += 0.5
        except KeyError: 
          pass # skip the first write when there is no data to write
        vars()[popName + "fastaSeq"] = list(fastaRef[CHR])

      # replace genotypes in fasta sequence if it differs from the REF
      if fastaRef[CHR][POS-1] != words[2]:
        sys.exit('Error: There is a discordance between fasta file and the REF column.\nFasta: %s\nTab: %s' % ([fastaRef.name[CHRnum-1], POS, fastaRef[CHR][POS-1]], words[0:3]))
      elif fastaRef[CHR][POS-1] != mostFreqGT and mostFreqGT != 'NA':
        vars()[popName + "fastaSeq"][POS-1] = mostFreqGT

    # to keep line tracking
    CHRprev = CHR
    fastaNum += int(fastaNumAdd)
    fastaNumAdd = 0
    
    # write any sequence that is not exists the tab file.
    while CHRnum > int(fastaRef.num[fastaNum]):
      fileoutput.write(">%s\n%s\n" % (fastaRef.name[fastaNum], fastaRef.sequence[fastaNum]))
      fastaNum += 1
      
    CHRnum = int(words[0].split("_")[1])

# write the last sequences that was modified
for popName in groups:
  fastaSeqP = ''.join(str(w) for w in vars()[popName + "fastaSeq"])
  fastaNameP = fastaRef.name[CHRnum-1]
  fileoutput.write(">%s.%s\n%s\n" % (fastaNameP, popName[-1], fastaSeqP))

# write any left sequence that were not modified by tab file:
while int(len(fastaRef.num)) > int(fastaRef.num[fastaNum]):
    fastaNum += 1
    fileoutput.write(">%s\n%s\n" % (fastaRef.name[fastaNum], fastaRef.sequence[fastaNum]))

fileoutput.close()
callsFile.close()

############################## to improve ##############################
"""
Code to use indels as well.
"""
