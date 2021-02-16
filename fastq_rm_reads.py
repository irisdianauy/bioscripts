#! /usr/bin/env python

"""
author @irisdianauy

This script:
[1] cleans up fastq files by removing reads based on read ID;
[2] can perfom this in parallelized batches;
[3] is faster than other methods but is memory intensive.
    Caution is needed when fastq files are large.
    This script is ideal to run in a large-memory HPC environment.

Future improvements
This currently assumes that the source and blacklist fastq files
have the same name.  In future edits, users will be allowed
to provide file patterns in a config file for cases where source
and blacklist fastq files have different names but follow a pattern.

Usage:
fastq_rm_reads.py </path/to/src/dir> </path/to/rmv/dir> </path/to/out/dir> <cores>
"""


import re
import subprocess
from os import path
from glob import glob
from sys import argv, exit
from multiprocessing import Pool


def stats_init(pstats):
    with open(pstats, "a+") as fstats:
        print("orig_fq_file,orig_reads,cleaned_fq_file,cleaned_reads",
              file=fstats)


def get_fastqs(pin, sfq):
    # All of the sfq files in pin
    lfqs = glob(f"{pin}/*.{sfq}")
    # All sample names present in lfqs
    lsamps = [(path.basename(p)).replace(sfq, "") for p in lfqs]
    # Dict with sample name to list of path
    dsamps = {samp: [s for s in lfqs if samp in s] for samp in lsamps}
    return(dsamps)


def get_fq_ids(pfq):
    lids = []
    if path.getsize(pfq) > 0:
        scmd = f"grep '^@' {pfq}"
        btraw = subprocess.check_output(scmd, shell=True)
        sids = str(btraw, "utf-8").replace("@", "")
        lids = sids.split("\n")
    return(lids)


def collate_fq_ids(lfqs):
    lids_all = []
    for pfq in lfqs:
        lids_all.extend(get_fq_ids(pfq))
    return(lids_all)


def fq_by_block(pfq):
    sfq_blk = "(?P<blk>^@(?P<id>[^\n]+)\n^[^\n]+\n^\+[^\n]*\n[^\n]+\n)"
    ore = re.compile(sfq_blk, re.M)
    with open(pfq, "r") as fsrc:
        ssrc = fsrc.read()
    ltblks = ore.findall(ssrc)
    return(ltblks)


def clean_fq(lids_rmv, psrc_fq, pclean):
    dblks = {t[1]: t[0] for t in fq_by_block(psrc_fq)}
    set_dblks = set(dblks.keys())
    set_diff_ids = set_dblks - set(lids_rmv)
    lblks = [dblks[sid] for sid in sorted(set_diff_ids)]

    with open(pclean, "w") as fclean:
        print("".join(lblks), file=fclean)

    with open(pstats, "a+") as fstats:
        psrc_fqb = path.basename(psrc_fq)
        pclean_b = path.basename(pclean)
        iblks = len(set_dblks)
        irmvs = len(lblks)
        print(f"{psrc_fqb},{iblks},{pclean_b},{irmvs}", file=fstats)


def clean_by_pair(lfq_pair):
    lsrcs, lrmvs = lfq_pair[0], lfq_pair[1]
    lids_rmv = collate_fq_ids(lrmvs)
    for psrc_fq in lsrcs:
        sbase = (path.basename(psrc_fq)).replace(".fq", "")
        pclean = f"{pout}/{sbase}.clean.fastq"
        clean_fq(lids_rmv, psrc_fq, pclean)


def clean_fqs(dsrcs, drmvs, pout, ip):
    lfq_pairs = [(dsrcs[s], drmvs[s]) for s in dsrcs]

    with Pool(ip) as p:
        p.map(clean_by_pair, lfq_pairs)


def main(psrc, prmv, pout, ip, pstats):
    stats_init(pstats)

    sfq = "fq"

    # Get source fastqs
    dsrcs = get_fastqs(psrc, sfq)

    # Get fastqs to remove
    drmvs = get_fastqs(prmv, sfq)

    # Produce cleaned fq/s
    clean_fqs(dsrcs, drmvs, pout, ip)


# """
if __name__ == "__main__":
    # Checkpoint
    try:
        psrc = path.abspath(argv[1])
        prmv = path.abspath(argv[2])
        pout = path.abspath(argv[3])
        ip = int(argv[4])
        pstats = f"{pout}/stats.csv"
    except IndexError:
        ppy = path.basename(__file__)
        shelp = f"\nUsage:\n{ppy} </path/to/src/dir> </path/to/rmv/dir> </path/to/out/dir> <cores>\n"
        print(shelp)
        exit()

    main(psrc, prmv, pout, ip, pstats)
# """
