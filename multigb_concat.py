#! /usr/bin/env python

"""
author @irisdianauy

The script converts a multi-genbank file to a single-record
genbank file (concatenated, numbering following the
gb record order).

Input: Multi-Genbank
"""


from os import path
from sys import argv, exit
from Bio import SeqIO, SeqFeature


def new_loc(istart, iend, istrand=None):
    ostr = SeqFeature.ExactPosition(istart)
    oend = SeqFeature.ExactPosition(iend)
    oloc = SeqFeature.FeatureLocation(ostr, oend, istrand)
    return(oloc)


def get_feature(istart, iend, sname, stype, istrand=None):
    ofeat_loc = new_loc(istart, iend, istrand)
    ofeat = SeqFeature.SeqFeature(ofeat_loc, type=stype)
    ofeat.qualifiers["locus_tag"] = sname
    return(ofeat)


def get_feature_contig(istart, iend, scont_name):
    ofeat = get_feature(istart, iend, scont_name, "contig")
    return(ofeat)


def main(pgbm, pgbc, sid):
    # Open the multi-genbank file
    ogbs = SeqIO.parse(pgbm, "genbank")

    # Initialize record for concatenated gb
    ogbc = next(ogbs)
    # Add contig feature
    ofeat_contig1 = get_feature_contig(0, len(ogbc), ogbc.name)
    ogbc.features.append(ofeat_contig1)
    # Change name
    ogbc.name = sid
    ogbc.id = sid
    ogbc.description = ogbc.description + ", concatenated"

    # Iterate over the rest of the genbank entries (start with 2)
    for ogb in ogbs:
        istr = len(ogbc)
        iend = istr + len(ogb)
        ogbc.seq += ogb.seq  # length of ogbc changes after this

        # Add contig feature
        ofeat_contig = get_feature_contig(istr, iend, ogb.name)
        ogbc.features.append(ofeat_contig)

        # Iterate over ogb features, edit, add to ogbc
        for ofeat in ogb.features:
            istr_new = ofeat.location.start.position + istr
            iend_new = ofeat.location.end.position + istr
            oloc_new = new_loc(istr_new, iend_new, ofeat.strand)
            ofeat.location = oloc_new
            ogbc.features.append(ofeat)

    # Write concatenated record to file
    SeqIO.write(ogbc, pgbc, "genbank")


# """
if __name__ == "__main__":
    try:
        pgbm = path.abspath(argv[1])
        pgbc = path.abspath(argv[2])
        sid = argv[3]
    except IndexError:
        ppy = path.basename(path.realpath(__file__))
        shelp = f"Usage\n{ppy} </path/to/multi.gbk> </path/to/concat.gbk> <id>"
        print(shelp)
        exit()

    main(pgbm, pgbc, sid)
# """
