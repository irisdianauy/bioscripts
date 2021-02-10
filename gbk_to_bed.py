#! /usr/bin/env python

"""
author @irisdianauy

The script converts a multi-genbank file into
a single bed file.

Feature: Specify features to extract (allows all).
"""


import pandas as pd
from os import path
from sys import argv, exit
from Bio import SeqIO


def get_id(ofeat):
    sid = (ofeat.qualifiers.get("db_xref", []))[0]
    sid = sid.replace("SEED:", "")
    return(sid)


def get_cds_info(ofeat):
    bpos = ofeat.strand > 0
    dcds = {"str": ofeat.location.start.position,
            "end": ofeat.location.end.position,
            "name": get_id(ofeat),
            "dxn": "+" if bpos else "-"}
    return(dcds)


def get_all_cds(ogb, lfeats):
    lcds_all = []
    for orec in ogb:
        dcds_ = {"chr": orec.id, "score": 1000}
        for ofeat in orec.features:
            if ofeat.type in lfeats:
                dcds = {**get_cds_info(ofeat), **dcds_}
                lcds_all.append(dcds)
    return(lcds_all)


def write_bed(lcds):
    odf = pd.DataFrame(lcds)
    odf.to_csv(pbed, sep="\t", index=False, header=False,
               columns=["chr", "str", "end", "name", "score", "dxn"])


def main(pgb, pbed, lfeats):
    # Open genbank file
    ogb = SeqIO.parse(pgb, "genbank")

    # Get all CDSs + info
    lcds_all = get_all_cds(ogb, lfeats)

    # Write locations in BED format
    write_bed(lcds_all)


# """
if __name__ == "__main__":
    try:
        pgb = path.abspath(argv[1])
        pbed = path.abspath(argv[2])
        lfeats = argv[3:]
    except IndexError:
        ppy = path.basename(path.abspath(__file__))
        shelp = f"Usage:\n{ppy} <input.gbk> <output.bed> <type1 type2 typen>"
        print(shelp)
        exit()

    main(pgb, pbed, lfeats)
# """
