# Dataset creation pipeline

Pipeline to create a dataset of closely-related bacterial isolates. The pipeline:

- downloads all RefSeq sequences from a given species.
- performs MLST with [tseemann/mlst](https://github.com/tseemann/mlst)
- and for a given desired sequence type and reference assembly:
  - evaluates mash distance of the chromosome to the reference assembly chromosome
  - selects all isolates with the specified sequence type, and at a mash distance smaller than a given threshold to the reference
  - exports all chromosomes from the dataset, along with metadata

## setup

The pipeline requires a working installation of conda/mamba and snakemake (the pipleine was designed for v8.10.8).
For cluster execution, also the `snakemake-executor-plugin-slurm` needs to be installed.
For convenience we provide the `snakemake_env.yml` from which the environment can be initialized.

## picking the dataset

The config file `config.yml` contains a list of the datasets to download:

```yml
species:
  saureus: "staphylococcus aureus"

mlst_scheme:
  saureus: "saureus"

dsets:
  saureus_ST5:
    species: "saureus"
    strain: "5"
    ref_acc: "GCF_000009665.1"
    threshold: 0.005
```

In this case `saureus_ST5` is a new dataset label, the species corresponds to _Staphylococcus aureus_, as indicated in the `species` entry.
The desired sequence type is ST5, as indicated by `strain`, with corresponding reference assembly `GCF_000009665.1`. All assemblies of this sequence type and at a mash distance of less than `0.005` will be retained in the dataset.

## running the pipeline

All of the datasets can be downloaded with:

```
snakemake all
```