#!/usr/bin/env python3
import numpy as np
import pandas as pd


df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
dflengths = pd.DataFrame(df.iloc[:, -1])
dflengths.index.name = "Chromosome"
dflengths.columns = ["E-P Distance"]
dflengths = dflengths[dflengths["E-P Distance"] >= 0]
dflengthsClean = dflengths[dflengths["E-P Distance"] <= 200000]

final = dflengthsClean.groupby("Chromosome")["E-P Distance"].mean()

final.to_csv(snakemake.output[0], sep="\t")
