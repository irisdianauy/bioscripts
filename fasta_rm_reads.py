#! /usr/bin/env python

"""
author @irisdianauy

This script removes contigs from a source fasta file,
based on contigs in one or many fasta files, collected in a single directory.
Requires that contigs to be removed have the same ID
as their counterpart in the source file.

Fast, even for very large fasta files, without using non-standard libraries.

This script can be useful in:
- creating a fasta file of contigs without identified contaminants
- creating a fasta file of contigs without known sequences
(e.g., Bacteria, Archaea, Eukarya, or Virus) so that the unknowns
can be explored further.
"""


from sys import argv, exit
from glob import glob
from os import path
import subprocess
import re


def stats_init(pstats):
    with open(pstats, "w") as fstats:
        print("orig_fa_file,orig_reads,cleaned_fa_file,cleaned_reads",
              file=fstats)


def get_fa_ids(pfa):
    lids = []
    if path.getsize(pfa) > 0:
        scmd = f"grep '^>' {pfa}"
        btraw = subprocess.check_output(scmd, shell=True)
        sids = str(btraw, "utf-8").replace(">", "")
        lids = sids.splitlines()
    return(lids)


def collate_fa_ids(lfas):
    lids_all = []
    for pfa in lfas:
        if path.exists(pfa):
            lids_all.extend(get_fa_ids(pfa))
    return(lids_all)


def fa_by_block(pfa):
    sfa_blk = "(?P<blk>^>(?P<id>[^\n]+)\n^[^>]+\n)"
    ore = re.compile(sfa_blk, re.M)
    with open(pfa, "r") as fsrc:
        ssrc = fsrc.read()
    ltblks = ore.findall(ssrc)
    return(ltblks)


def clean_fa(lids_rmv, psrc, pout):
    dblks = {t[1]: t[0] for t in fa_by_block(psrc)}
    set_dblks = set(dblks.keys())
    set_diff_ids = set_dblks - set(lids_rmv)
    lblks = [dblks[sid] for sid in sorted(set_diff_ids)]

    with open(pout, "w") as fout:
        print("".join(lblks), file=fout)

    with open(pstats, "a+") as fstats:
        psrcb = path.basename(psrc)
        pout_b = path.basename(pout)
        iblks = len(set_dblks)
        irmvs = len(lblks)
        print(f"{psrcb},{iblks},{pout_b},{irmvs}", file=fstats)


def main(psrc, prmv, pout, pstats):
    stats_init(pstats)

    # Get fastas to remove
    lrmvs = glob(f"{prmv}/*.fasta")

    # Get IDs to remove
    lids_rmv = collate_fa_ids(lrmvs)

    # Produce cleaned fa/s
    clean_fa(lids_rmv, psrc, pout)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        psrc = path.abspath(argv[1])
        prmv = path.abspath(argv[2])
        pout = path.abspath(argv[3])
        pstats = f"{pout}.stats"
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"\nUsage:\n{ppy} </path/to/src/fasta> </path/to/rmv/dir> </path/to/out/fasta>\n"
        print(shelp)
        exit()

    main(psrc, prmv, pout, pstats)
# """
