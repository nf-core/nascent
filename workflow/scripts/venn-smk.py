"""
Takes in 3 files
0 - the left
1 - the right
2 - the overlap
"""
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import os

with open(snakemake.input[0]) as f:
    count1 = sum(1 for _ in f)

with open(snakemake.input[1]) as f:
    count2 = sum(1 for _ in f)

with open(snakemake.input[2]) as f:
    common = sum(1 for _ in f)

fileName0 = os.path.basename(os.path.normpath(snakemake.input[0]))
fileName1 = os.path.basename(os.path.normpath(snakemake.input[1]))
fileName2 = os.path.basename(os.path.normpath(snakemake.input[2]))

v = venn2(
    subsets=((count1 - common), (count2 - common), common),
    set_labels=(fileName0, fileName1, fileName2),
)
v.get_patch_by_id("100").set_color("blue")
v.get_patch_by_id("010").set_color("red")
v.get_patch_by_id("110").set_color("purple")

plt.title(snakemake.params["title"])

plt.savefig(snakemake.output[0])