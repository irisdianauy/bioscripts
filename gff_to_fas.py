#! /usr/bin/env python

"""
author @irisdianauy

This script reads a GFF file
and a separate scaffolds fasta file,
and extracts all mRNAs.

The script is suitable for use on GFF files
that do not contain sequences (i.e., no ##FASTA lines),
but separate fasta files are available.

Output is a fasta file of all predicted mRNAs.
"""


import re
import pandas as pd
from os import path
from sys import argv, exit
from Bio import SeqIO


def get_feature_id(sattr):
    sid = re.match("ID=([^;]+);", sattr).groups()[0]
    return(sid)


def gff_to_df(pgff):
    """
    Generate a dataframe from a GFF file.
    """
    lcols_use = [0, 2, 3, 4, 8]
    lcols_names = ["seqname", "feature", "start", "end", "attribute"]
    odf = pd.read_csv(pgff, sep="\t", header=None, index_col=None,
                      skiprows=[0], usecols=lcols_use,
                      names=lcols_names)
    odf["fid"] = odf.attribute.apply(get_feature_id)
    odf.set_index("fid", inplace=True)
    return(odf)


def extract_seq(trow, dscaffs):
    sid = trow[0]
    scaff_id = trow[1].seqname
    istart = int(trow[1].start)-1
    iend = int(trow[1].end)
    oseq_scaff = dscaffs[scaff_id]
    oseq = oseq_scaff[istart:iend]
    oseq.id = sid
    oseq.name = ""
    oseq.description = ""
    return(oseq)


def write_feature_seqs(odf, dscaffs, pfas):
    # In these gffs, there are only 'mRNA' and 'CDS'
    sfeature = "mRNA"
    odf_feat = odf.loc[odf["feature"] == sfeature, :]
    with open(pfas, "w") as fas:
        for trow in odf_feat.iterrows():
            oseq = extract_seq(trow, dscaffs)
            SeqIO.write(oseq, fas, "fasta")


def main(pgff, pscaffs, pfas):
    odf = gff_to_df(pgff)
    dscaffs = SeqIO.to_dict(SeqIO.parse(pscaffs, "fasta"))
    write_feature_seqs(odf, dscaffs, pfas)


# """
if __name__ == "__main__":
    try:
        pgff = path.abspath(argv[1])
        pscaffs = path.abspath(argv[2])
        pfas = path.abspath(argv[3])
    except IndexError:
        ppy = path.basename(path.realpath(__file__))
        shelp = f"Usage:\n{ppy} <input_gff> <input_scaffolds> <output_fasta>"
        print(shelp)
        exit()

    main(pgff, pscaffs, pfas)
# """
