# gene_search_in_encode
Examine transcript variety of gene of interest in ENCODE long-read RNA sequencing data


## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/HYzhang800/gene_search_in_encode.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd gene_search_in_encode
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate encode_gene_search`.

It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
cd gene_search_in_encode
conda activate encode_gene_search
snakemake --use-conda -n all
```
Run the pipeline using the following command:

```bash
snakemake --use-conda -j <num_cores> all
```


You can visualise the processes to be executed in a DAG:

```bash
snakemake --dag | dot -Tpng > dag.png
```

To exit a running `snakemake` pipeline, hit `ctrl+c` on the terminal. If the pipeline is running in the background, you can send a `TERM` signal which will stop the scheduling of new jobs and wait for all running jobs to be finished.

```bash
killall -TERM snakemake
```

To deactivate the conda environment:
```bash
conda deactivate
```
