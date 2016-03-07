#!/usr/bin/env python

"""
Created: 2016-03-05

@author: irisdianauy

This tool extractS CDS position info from NCBI features table/s.
The output is an exon file for use in trainGlimmerHMM.

Future edits:
[ 1] Provide support for a multi-Genbank file.
[ 2] Provide support for hybrid input: feat table and multi-Genbank
[ 3] Provide option to specify inclusion of fuzzy starts/ends.
"""

from __future__ import print_function
import os.path
import argparse
import re


#"""
# Define regular expressions
# Regex that matches a whole feature block in an NCBI features table
reFeatBlock = re.compile(">Feature[\s][\d\|\w\.]+\n(?:(?:[\t]|[\>\<](?=\d)|[\d])[ \t\S]+\n)+")

# Regex for capturing the ID of a whole feature block
reId = re.compile(">Feature[\s](?P<id>[\S]+)\n")

# Regex that matches a CDS block in a feature block
reCdsBlock = re.compile("^[\d]+\t[\d]+\tCDS\n(?:[\d]+\t[\d]+\n)+(?=\t)", re.M)

# Regex that matches a CDS block in a feature block, with fuzz
reCdsBlkfz = re.compile("^[\>\<\d]+\t[\>\<\d]+\tCDS\n(?:[\>\<\d]+\t[\>\<\d]+\n)+(?=\t)", re.M)

# Regex that matches an exon pair in a CDS block
rePosPair = re.compile("(?P<start>[\S]+)\t(?P<stop>[\S]+)")
#"""


#"""
# Define classes & functions
def getExons(sId, lCoords):
    """This function assumes that lCoords is a list of dictionaries
    containing the positions of exons in a gene."""
    sGene = ""
    for dCoords in lCoords:
        dCoords["id"] = sId
        sGene += "{id} {start} {stop}\n".format(**dCoords)
    return sGene

class oFeatBlock:
    def __init__(self, sLines):
        rsId = reId.search(sLines)
        dId = rsId.groupdict()
        self.id = dId["id"]
        self.lfeats = reCdsBlock.findall(sLines)
        self.numfeats = len(self.lfeats)
    def printFeats(self, fExon):
        for sFeat in self.lfeats:
            lExons = [e.groupdict() for e in rePosPair.finditer(sFeat)]
            print(getExons(self.id, lExons), file=fExon)
#"""


#"""
# Set and parse command line arguments
sDescTool = '''This tool extracts CDS positions from 
one or more Genbank feature tables and generates
an exon file suitable for trainGlimmerHMM.'''

sHelpIn = '''Input Genbank features table file/s.
For multiple files, specify this option multiple times.
This option is required.'''

sHelpOut = '''Name of new file in which results are written.
Defaults to the name of the input file suffixed with '.exon'.
This option is required.'''

oOpts = argparse.ArgumentParser(description = sDescTool)
oOpts.add_argument("-i", "--infeat",
                   metavar = "file.features",
                   action = "append",
                   help = sHelpIn)
oOpts.add_argument("-o", "--outexon",
                   metavar = "file.exons",
                   help = sHelpOut)
oParsedOpts = oOpts.parse_args()
#"""


#"""
# Checkpoints
# Check if options which can't be set as positional args are supplied
# Check supplied options for conflict
print("\nChecking options now...")
dReqd = {"sInfeat" : [oParsedOpts.infeat, "input features table/s"],
         "sOutexon": [oParsedOpts.outexon, "output exons file path/name"]}

def checkOpts(lParam):
    if lParam[0]:
        if isinstance(lParam[0], list):
            return len(lParam[0])
        else:
            return 1
    else:
        print("No {} provided.".format(lParam[1]))
        return 0

for sReqd in dReqd:
    iReqd = checkOpts(dReqd[sReqd])
    if iReqd == 0:
        print("Supply all required inputs.")
        break
    else:
        dReqd[sReqd].append(iReqd)

# Check if outfile exists
bOutExists = os.path.isfile(oParsedOpts.outexon)
if bOutExists:
    print("Output filename exits. Supply a new one.")

# Include all conditions in a compound checkpt for single check
bProceed = False
if iReqd != 0 and not bOutExists:
    bProceed = True
#"""


#"""
# Proceed if all things check out
if bProceed:
    print("Requirements supplied. Proceeding.")
    with open(oParsedOpts.outexon, "w") as fOutExon:
        for pIn in oParsedOpts.infeat:
            fIn = open(pIn).read()
            sIn = str(fIn)
            lFeats = reFeatBlock.findall(sIn)
            for sFeats in lFeats:
                oFeat = oFeatBlock(sFeats)
                oFeat.printFeats(fOutExon)
else:
    pass
#"""


# ==============================================================================
#"""
# Regular expressions: Notes
# Regex that matches a whole feature block:
# ">Feature[\s][\d\|\w\.]+\n(?:(?:[\t]|[\>\<](?=\d)|[\d])[ \t\S]+\n)+"
## Broken down into parts:
#### First line of the block
#### ">Feature[\s][\d\|\w\.]+\n"
#### Succeeding lines
#### "(?:(?:[\t]|[\>\<](?=\d)|[\d])[ \t\S]+\n)+"
###### where 
###### "[\t]|[\>\<](?=\d)|[\d])"
###### specifies the start of every line that isn't the first line, and
###### "[ \t\S]+"
###### is the pattern for all succeeding characters
