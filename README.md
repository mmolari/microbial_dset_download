# Microbial dataset download pipeline

Pipeline to download and characterize complete bacterial genomes from NCBI RefSeq.

Given a series of species, to be specified in the `config.yml` file, the pipeline:

- downloads all `complete` RefSeq sequences from a given species.
- performs MLST with [tseemann/mlst](https://github.com/tseemann/mlst)
- builds an approximate phylogeny using neighbour joining from mash distances using [attotree](https://github.com/karel-brinda/attotree)
- produces a summary metadata table and figures.
- optionally, by running the `cluster_all` target, clusters genomes creating clades of closely related isolates with two methods:
  - selecting internal nodes of the approximate phylogenetic tree that are associated to a particular Sequence Type.
  - clustering the genomes using HDBSCAN on the mash distances, picking a threshold parameter that makes the clustering compatible with the MLST scheme used.
- optionally, running the `export_all` target creates an `export` directory with:
  - for each species, a directory with the compressed genomes, the metadata table, the phylogenetic tree, MLST results, compressed mash distances, anb the HDBSCAN clustering results.
  - and a global `all_clusters.yaml` file with a list of all cluster names per species.

## setup

The pipeline requires a working installation of conda/mamba and snakemake (the pipeline was designed for v8.10.8).
For cluster execution, also the `snakemake-executor-plugin-slurm` needs to be installed.
For convenience we provide the `snakemake_env.yml` from which the environment can be initialized.

Optionally you can add your ncbi api key in a file named `ncbi_api_key.txt`. This is not required, but will speed up the download process. You can obtain an API key from [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/).

## picking the dataset

The config file `config.yml` contains a list of the datasets to download:

```yml
species:
  saureus: "staphylococcus aureus"

mlst_scheme:
  saureus: "saureus"
```

In this case the pipeline will download the `saureus` dataset, corresponding to the _Staphylococcus aureus_ species.
In the `species` dictionary you can specify for each dataset the search term for the taxonomy search in NCBI.
The `mlst_scheme` dictionary specifies the MLST scheme to use for each dataset, which
should correspond to one of the [available schemes](https://github.com/tseemann/mlst?tab=readme-ov-file#available-schemes) in the `tseemann/mlst` repository.

## running the pipeline

All of the datasets can be downloaded with:

```sh
snakemake all
```

Or for SLURM cluster execution:

```sh
snakemake all --profile profiles/slurm
```

Substitute `all` with `cluster_all` to run the clustering step, which will create a directory `clusters/` with the clustering results, or with `export_all` to create the `export` directory with all the results.