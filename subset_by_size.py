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


def gen_new_file(pfile, pnew, imin, sformat):
    orecs = SeqIO.parse(pfile, sformat)
    lpass = [orec for orec in orecs if len(orec) >= imin]
    with open(pnew, "w") as fnew:
        for orec in lpass:
            SeqIO.write(orec, fnew, sformat)


def main(pin, pout, imin, sformat):
    lfiles = glob(f"{pin}/*.*")
    for pfile in lfiles:
        sfile, sext = path.splitext(pfile)
        samp = path.basename(sfile)
        pnew = f"{pout}/{samp}_{imin}up{sext}"
        gen_new_file(pfile, pnew, imin, sformat)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        pin = path.abspath(argv[1])
        pout = path.abspath(argv[2])
        imin = int(argv[3])
        sformat = argv[4]
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"Usage:\t{ppy} </path/to/in/fasta/dir> </path/to/outdir/with/prefix> <min_len> <fasta|genbank>"
        print(shelp)
        exit()

    main(pin, pout, imin, sformat)
# """
