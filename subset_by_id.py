#! /usr/bin/env python

"""
author @irisdianauy

This script subsets a fasta file
using a whitelist of IDs.

"""


import re
from sys import argv, exit
from os import path


def fa_by_block(pfa):
    """
    From a multifasta file,
    return a list of fasta string blocks.
    """
    sfa_blk = "(?P<blk>^>(?P<id>[^\n]+)\n^[^>]+\n)"
    ore = re.compile(sfa_blk, re.M)
    with open(pfa, "r") as fsrc:
        ssrc = fsrc.read()
    ltblks = ore.findall(ssrc)
    return(ltblks)


def clean_fa(lids_sub, psrc, pout):
    """
    From a source multifasta file,
    select fasta blocks by ID
    and write selection to a new file.
    """
    dblks = {t[1]: t[0] for t in fa_by_block(psrc)}
    lblks = [dblks[sid] for sid in lids_sub]
    with open(pout, "w") as fout:
        print("".join(lblks), file=fout)


def main(psrc, psub, pout, pstats):
    # Generate a list of ids
    lids_sub = open(psub).read().splitlines()

    # Produce subset fasta
    clean_fa(lids_sub, psrc, pout)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        psrc = path.abspath(argv[1])
        psub = path.abspath(argv[2])
        pout = path.abspath(argv[3])
        pstats = f"{pout}/stats.csv"
    except IndexError:
        shelp = "\nUsage:\nsubset_by_id.py </path/to/src/fasta> </path/to/subset/ids> </path/to/out/fasta>\n"
        print(shelp)
        exit()

    main(psrc, psub, pout, pstats)
# """
