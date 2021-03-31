#!/usr/bin/env python
import pandas as pd
import argparse


def bed2saf(bed, output):
    bed = pd.read_csv(bed, sep="\t")
    bed.columns = ["Chr", "Start", "End", "GeneID", "frame", "Strand"]
    # TODO Clever naming scheme
    bed["GeneID"] = (
        bed["Chr"] + "-" + bed["Start"].astype(str) + "-" + bed["End"].astype(str)
    )
    bed.to_csv(
        output,
        sep="\t",
        columns=["GeneID", "Chr", "Start", "End", "Strand"],
        index=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Convert a custom bed
        to a SAF annotation."""
    )
    parser.add_argument("bed", type=str, help="Custom transgene sequence")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="transgenes.saf",
        type=str,
        help="SAF output",
    )
    args = parser.parse_args()
    bed2saf(args.bed, args.output)
