"""
Takes in 3 files
0 - the left
1 - the right
2 - the overlap
"""

from matplotlib_venn import venn2

with open(snakemake.input[0]) as f:
    count1 = sum(1 for _ in f)

with open(snakemake.input[1]) as f:
    count2 = sum(1 for _ in f)

with open(snakemake.input[2]) as f:
    common = sum(1 for _ in f)


venn2(
    subsets=((count1 - common), (count2 - common), common),
    set_labels=(snakemake.input[0], snakemake.input[1], snakemake.input[2]),
)
v.get_patch_by_id("100").set_color("blue")
v.get_patch_by_id("010").set_color("red")
v.get_patch_by_id("110").set_color("purple")

plt.title(snakemake.params["title"])

plt.savefig(snakemake.output[0])
