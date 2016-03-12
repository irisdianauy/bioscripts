#!/usr/bin/env python
# This script is for summarizing BLAST XML results into a table

# Edits
# [ 1] Changed from comma to tab as delimiter file
# [ 2] Lengthened the hit definition line
# [ 3] Included check of BLAST program used to control calculation of qcov
# Future edit:
# [ 1] DONE: Consider limits to hits and hsps to display
# [ 2] DONE: Determine BLAST program used and use this to control calculation on query coverage
# [ 3] Rewrite to make it better, more readable; use functions.


#"""
import argparse
from Bio.Blast import NCBIXML
import re
#"""


#"""
options = argparse.ArgumentParser(description="Create a summary from BLAST results")
options.add_argument("input",
                     metavar="blast_result.xml",
                     help="BLASTX result in XML format")
options.add_argument("-o", "--output",
                     metavar="blast_result.sum",
                     help="Output filename. Extension '.sum' is appended by default. Defaults to <input>.sum")
options.add_argument("-n", "--nohits",
                     metavar="blast_result.nohits",
                     help="Output filename. Extension '.nohits' is appended by default. Defaults to <input>.nohits, if no filename is specified")
options.add_argument("-t", "--max_hits",
                     metavar='int',
                     type=int,
                     help="Integer. Maximum hits per query to be displayed. Defaults to all hits.")
options.add_argument("-p", "--max_hsps",
                     metavar='int', type=int,
                     help="Integer. Maximum HSPs per query hit to be displayed. Defualts to all HSPs.")
parsedOpts = options.parse_args()
print options
print parsedOpts
#"""


#"""
# For input file
blast = NCBIXML.parse(open(parsedOpts.input), "rU")

reProg = re.compile(
    "<BlastOutput_program>(?P<prog>[\w]*blast[\w]*)</BlastOutput_program>")
with open(parsedOpts.input) as fBlast:
    for sLine in fBlast:
        if reProg.search(sLine):
            rsProg = reProg.search(sLine)
            dProg = rsProg.groupdict()
            sProg = dProg["prog"]
            break
        else:
            pass

# For output file
if parsedOpts.output:
    if (parsedOpts.output).endswith(".sum"):
        outFileName = parsedOpts.output
    else:
        outFileName = "{}.sum".format(parsedOpts.output)
else:
    outFileName = "{}.sum".format((parsedOpts.input).split(".")[0])
outFile = open(outFileName, "w")
print >> outFile, "QUERY\tQUERY_LEN\tHIT_ACC\tHIT_SP\tHIT_DSC\tHSP_LEN\tHSP_EVAL\tHSP_QCOV\tHSP_ID\tHSP_SIM"

# For no-hits file
if parsedOpts.nohits:
    if (parsedOpts.nohits).endswith(".nohits"):
        pNohits = parsedOpts.nohits
    else:
        pNohits = "{}.nohits".format(parsedOpts.nohits)
else:
    pNohits = "{}.nohits".format((parsedOpts.input).split(".")[0])
fNohits = open(pNohits, "w")
#"""


#"""
for record in blast:
    query = record.query
    queryLen = record.query_length

    if record.alignments:

        currHitCount = 0

        # Record-related attributes
        for recHit in record.alignments:
            # Hit/Alignment-related attributes
            recHitAcc = recHit.accession
            if "[" in recHit.title:
                recHitSp = (recHit.title).split(
                    "[")[1].split("]")[0].replace(",", ";")
            else:
                recHitSp = "None stated in definition"
            # for recHitDsc
            recHitDef = (recHit.hit_def).replace(",", ";")
            if recHitDef.startswith("gi|"):
                recHitDef_ = (recHitDef.split(">")[0]).split("|")[-1]
                # Threshold in the following used to be 50
                if len(recHitDef_) <= 200:
                    recHitDsc = recHitDef_
                else:
                    recHitDsc = recHitDef_[:197] + "..."
            else:
                if len(recHitDef) <= 200:
                    recHitDsc = recHitDef
                else:
                    recHitDsc = recHitDef[:197] + "..."

            currHspCount = 0

            # HSP-related attributes
            for hitHsp in recHit.hsps:
                hitHspLen = hitHsp.align_length
                hitHspEval = hitHsp.expect
                if sProg == "blastx":
                    hitHspQcov = (
                        float(hitHsp.query_end - hitHsp.query_start + 1) / 3) / queryLen
                else:
                    # For blastn or blastp
                    hitHspQcov = float(
                        hitHsp.query_end - hitHsp.query_start + 1) / queryLen
                hitHspId = float(hitHsp.identities) / hitHspLen
                hitHspSim = float(hitHsp.positives) / hitHspLen

                # Printing to file
                lsMetrics = (query, queryLen, recHitAcc, recHitSp, recHitDsc,
                             hitHspLen, hitHspEval, hitHspQcov, hitHspId, hitHspSim)
                print >> outFile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:0.4f}\t{:0.4f}\t{:0.4f}".format(
                    *lsMetrics)

                currHspCount += 1
                if parsedOpts.max_hsps:
                    if currHspCount >= parsedOpts.max_hsps:
                        break

            currHitCount += 1
            if parsedOpts.max_hits:
                if currHitCount >= parsedOpts.max_hits:
                    break

    else:
        print >> fNohits, query
#"""

blast.close()
outFile.close()
fNohits.close()
