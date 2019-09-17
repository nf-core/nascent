import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in count data

df = pd.read_csv(snakemake.input[0], sep="\t", index_col="Geneid")
print(df)

# Get on the raw counts
dfCounts = df.filter(regex="\d", axis=1)

# Get the lengths of each gene
dfLength = df["Length"]


# Calculate reads per kilobase by taking
# ```
# counts / length
# ```
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
rpk = dfCounts.divide(dfLength, axis="index")

# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
rpkmFactor = rpk.sum(axis=0) / 1e6

# Divide the RPK values by the “per million” scaling factor. This gives you TPM.
dftpm = rpk / rpkmFactor


# Filter out all the zeros
dftpm = dftpm[dftpm > 0].dropna()

print(dftpm.head)

dftpm.to_csv(snakemake.output[0], sep="\t")
