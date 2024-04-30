#!/usr/bin/env python3

import platform
import pandas as pd

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

df_subread = pd.read_csv("$subread", sep="\\t", comment="#", index_col=0)
df_subread["Name"] = df_subread.index
for column in df_subread.columns[5:]:
    df_bed = df_subread[["Chr", "Start", "End", "Name", column, "Strand"]]
    sample_id = column[:-len(".sorted.bam")]
    df_bed.to_csv(f"{sample_id}.bed", sep="\\t", header=False, index=False)

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))