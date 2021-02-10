#! /usr/bin/env python

"""
Subset multifasta files by excluding sequences
that are less than a given minimum.

Works on multiple multifasta files in a dir.

Currently has Bio.SeqIo dependency.
"""


from os import path
from glob import glob
from sys import argv, exit
from Bio import SeqIO


def gen_new_file(pfas, pnew, imin):
    orecs = SeqIO.parse(pfas, "fasta")
    lpass = [orec for orec in orecs if len(orec) >= imin]
    with open(pnew, "w") as fnew:
        for orec in lpass:
            SeqIO.write(orec, fnew, "fasta")


def main(pin, pout, imin):
    lfas = glob(f"{pin}/*.fasta")
    for pfas in lfas:
        samp = (path.basename(pfas)).replace(".fasta", "")
        pnew = f"{pout}/{samp}_{imin}up.fasta"
        gen_new_file(pfas, pnew, imin)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        pin = path.abspath(argv[1])
        pout = path.abspath(argv[2])
        imin = int(argv[3])
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"Usage:\t{ppy} </path/to/in/fasta/dir> </path/to/outdir/with/prefix> <min_len>"
        print(shelp)
        exit()

    main(pin, pout, imin)
# """
