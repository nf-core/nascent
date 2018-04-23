# IMR90
## Getting Started
```
# clone workflow into working directory
git clone https://bitbucket.org/user/myworkflow.git path/to/workdir
cd path/to/workdir

# edit config and workflow as needed
vim config.yaml

# install dependencies into isolated environment
conda create -n myworkflow --file requirements.txt

# activate environment
source activate myworkflow

# execute workflow
snakemake -n
```
