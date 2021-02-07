#! /usr/bin/env python

"""
Split fasta file into files, <int> max per file.
Originally made for web-based eggnog-mapper submissions.
"""

from sys import argv, exit
from os import path
from Bio import SeqIO
from math import ceil


def popslice(lrecs, i):
    lcurr = lrecs[:i]
    del(lrecs[:i])
    return(lcurr)


def write_split_faa(pnew, lrecs, i):
    lcurr = popslice(lrecs, i)
    with open(pnew, "w") as fnew:
        for o in lcurr:
            SeqIO.write(o, fnew, "fasta")


def main(pin, pout, i):
    orecs = SeqIO.parse(pin, "fasta")
    lrecs = list(orecs)
    ilen = len(lrecs)
    if ilen <= i:
        print("No need to split.")
        exit()
    else:
        k = len(str(ceil(ilen/i)))
    sfile, sext = path.splitext(pin)
    samp = path.basename(sfile)
    j = 1
    while lrecs:
        sj = f"sub{j:0{k}}"
        pnew = f"{pout}/{samp}_{sj}{sext}"
        write_split_faa(pnew, lrecs, i)
        j += 1


if __name__ == "__main__":
    # Checkpoint
    try:
        pin = path.abspath(argv[1])
        pout = path.abspath(argv[2])
        i = int(argv[3])
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"Usage:\t{ppy} <input_fasta> <output_dir> <max_records>"
        print(shelp)
        exit()

    main(pin, pout, i)
