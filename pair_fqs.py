#! /usr/bin/env python

"""
author @irisdianauy

Pair unbalanced paired-end reads
(e.g., for products of QC).
"""


import re
from sys import argv, exit
from os import path
from collections import OrderedDict


# Set fq pattern
slines = "(?P<id>^@[^\s]+)[^\n]+\n"   # fq line 1
slines += "^[^@\n]+\n"                # fq line 2
slines += "^\+\n"                     # fq line 3
slines += "^[^\n]+\n"                 # fq line 4
sfq_blk = f"(?P<blk>{slines})"        # fq block
ore = re.compile(sfq_blk, re.M)


def fq_by_block(pfq):
    with open(pfq, "r") as fq:
        sfq = fq.read()
    ltblks = ore.findall(sfq)
    return(ltblks)


def get_unpaired_read_ids(pfq1, pfq2):
    ltblks1 = fq_by_block(pfq1)
    ltblks2 = fq_by_block(pfq2)
    set1 = set([t[1] for t in ltblks1])
    set2 = set([t[1] for t in ltblks2])
    set_rmv = set1 ^ set2
    return(ltblks1, ltblks2, set_rmv)


def write_block(lblks, pout_fq):
    with open(pout_fq, "w") as fout_fq:
        print("".join(lblks), file=fout_fq)


def pair_fq(ltblks, set_rmv, paired):
    dblks = OrderedDict({t[1]: t[0] for t in ltblks})
    lids = list(dblks.keys())
    lblks = [dblks[sid] for sid in lids
             if sid not in set_rmv]
    lblks_rmv = [dblks[sid] for sid in lids
                 if sid in set_rmv]
    write_block(lblks, paired)
    return(lblks_rmv)


def pair_fqs(ltblks1, ltblks2, set_rmv, samp, pout):
    lblks_rmv = []
    for i, ltblks in enumerate([ltblks1, ltblks2]):
        paired = f"{pout}/{samp}_R{i+1}.paired.fq"
        lblks_rmv += pair_fq(ltblks, set_rmv, paired)
    punpaired = f"{pout}/{samp}.unpaired.fq"
    write_block(lblks_rmv, punpaired)


def main(pfq1, pfq2, samp, pout):
    # Get IDs of unpaired reads
    ltblks1, ltblks2, set_rmv = get_unpaired_read_ids(pfq1, pfq2)

    # Produce paired read1 & read2 fq, and an unparied fq
    pair_fqs(ltblks1, ltblks2, set_rmv, samp, pout)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        pfq1 = path.abspath(argv[1])
        pfq2 = path.abspath(argv[2])
        samp = argv[3]
        pout = path.abspath(argv[4])
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"\nUsage:\n{ppy} </path/to/fq1> </path/to/fq2> <samp> </path/to/out/dir>\n"
        print(shelp)
        exit()

    main(pfq1, pfq2, samp, pout)
# """
