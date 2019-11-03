#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# Read in tpm values

df = pd.read_csv(snakemake.input[0], sep="\t", index_col="Geneid")

# Take log2 of samples
dflog = np.log2(df)

# Calculate percent change
dfFold = dflog.pct_change(axis="columns")
# Drop 0h since it won't have any fold change
dfFold.dropna(axis=1, inplace=True)

# Drop anything > 24h
try:
    dfFold = dfFold.drop(columns=["IMR48h", "IMR72h", "GM48h", "GM72h"])
except KeyError:
    pass  # do nothing!
# Drop anything that don't have a 1 fold change up or down(easy to change), might even make it a parameter in snakemake
dfDGEup = dfFold[dfFold >= snakemake.params["cutoff"]].dropna(thresh=1)
dfDGEdown = dfFold[dfFold <= -(snakemake.params["cutoff"])].dropna(thresh=1)
frames = [dfDGEup, dfDGEdown]
dfDGE = pd.concat(frames)

# This brings the TPM values back into the list of DE genes
dfDGE = df.merge(dfDGE, on="Geneid", suffixes=("", "_y"))
dfDGE = dfDGE.drop(dfDGE.filter(like="_y", axis=1).columns, axis=1)

dfDGE.to_csv(snakemake.output[0], sep="\t")
