#!/usr/bin/env python
import pandas as pd

bed = pd.read_csv(snakemake.input[0], sep='\t')
bed.columns = ["Chr", "Start", "End", "GeneID", "frame", "Strand"]
bed.to_csv(snakemake.output[0], sep='\t',columns=["GeneID", "Chr", "Start", "End", "Strand"],index=False)
