configfile: "config.yml"


localrules:
    download_ref,
    download_dataset,


ncbi_api_key = ""
try:
    with open("ncbi_api_key.txt", "r") as f:
        ncbi_api_key = f.read().strip()
except:
    print("No NCBI API key found. Optionally add your key in ncbi_api_key.txt")


species = config["species"]


rule download_ref:
    output:
        "data/ref/record/{ref_acc}.fa",
    conda:
        "envs/ncbi.yml"
    params:
        api=f"--api-key {ncbi_api_key}" if ncbi_api_key else "",
    shell:
        """
        datasets download genome accession {wildcards.ref_acc} \
            {params.api} \
            --filename {wildcards.ref_acc}.zip
        unzip {wildcards.ref_acc}.zip -d {wildcards.ref_acc}
        mv {wildcards.ref_acc}/ncbi_dataset/data/*/*.fna {output}
        rm -r {wildcards.ref_acc} {wildcards.ref_acc}.zip
        """


rule ref_to_chromosome:
    input:
        rules.download_ref.output,
    output:
        "data/ref/chromosome/{ref_acc}.fa",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/extract_chromosome.py --fa {input} --out_file {output}
        """


rule download_dataset:
    output:
        "data/species/{species}/ncbi.zip",
    params:
        ncbi_sp=lambda w: species[w.species],
        api=f"--api-key {ncbi_api_key}" if ncbi_api_key else "",
    conda:
        "envs/ncbi.yml"
    shell:
        """
        datasets download genome taxon '{params.ncbi_sp}' \
            {params.api} \
            --assembly-level complete \
            --assembly-source 'RefSeq' \
            --exclude-atypical \
            --mag exclude \
            --filename {output}
        """


rule expand_dataset:
    input:
        rules.download_dataset.output,
    output:
        ds=directory("data/species/{species}/ncbi"),
        info="data/species/{species}/info.jsonl",
    conda:
        "envs/ncbi.yml"
    shell:
        """
        TMP=$(mktemp -d)
        unzip {input} -d $TMP
        mkdir -p {output.ds}
        for f in $TMP/ncbi_dataset/data/*/*.fna; do
            mv $f {output.ds}/$(basename $(dirname $f)).fa
        done
        mv $TMP/ncbi_dataset/data/*.jsonl {output.info}
        rm -r $TMP
        """


rule info_to_tsv:
    input:
        rules.expand_dataset.output.info,
    output:
        "data/species/{species}/info.tsv",
    params:
        fields=",".join(config["metadata"]),
    conda:
        "envs/ncbi.yml"
    shell:
        """
        dataformat tsv genome \
            --fields {params.fields} \
            --inputfile {input} > {output}
        """


rule chromosome_fa:
    input:
        rules.expand_dataset.output.ds,
    output:
        directory("data/species/{species}/chromosomes"),
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        mkdir -p {output}
        for f in {input}/*.fa; do
            python3 scripts/extract_chromosome.py --fa $f --out_dir {output}
        done
        """


rule assembly_to_chrom_acc:
    input:
        rules.expand_dataset.output.ds,
    output:
        "data/species/{species}/assembly_to_chrom.tsv",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_to_chrom_acc.py --fas {input} --out {output}
        """


rule attotree:
    input:
        fas=rules.chromosome_fa.output,
    output:
        tree="data/species/{species}/attotree.nwk",
    conda:
        "envs/attotree.yml"
    threads: 4
    params:
        k=21,
        s=50000,
    shell:
        """
        attotree {input.fas}/*.fa \
            -o {output.tree} \
            -t {threads} \
            -k {params.k} \
            -s {params.s}
        """


rule mlst:
    input:
        rules.chromosome_fa.output,
    output:
        "data/species/{species}/mlst.tsv",
    params:
        scheme=lambda w: config["mlst_scheme"][w.species],
    conda:
        "envs/mlst.yml"
    shell:
        """
        mlst \
            --scheme {params.scheme} \
            {input}/*.fa > {output}
        """


# rule fig_mash_dist:
#     input:
#         mash=lambda w: expand(
#             rules.mash_dist.output,
#             species=dsets[w.dset]["species"],
#             ref_acc=dsets[w.dset]["ref_acc"],
#         ),
#         mlst=lambda w: expand(
#             rules.mlst.output,
#             species=dsets[w.dset]["species"],
#         ),
#     output:
#         "results/datasets/{dset}/figs/mash_dist.pdf",
#     params:
#         st_num=lambda w: dsets[w.dset]["strain"],
#         thr=lambda w: dsets[w.dset]["threshold"],
#     conda:
#         "envs/bioinfo.yml"
#     shell:
#         """
#         python3 scripts/plot_mash_dist.py \
#             --mash {input.mash} \
#             --mlst {input.mlst} \
#             --ST_num {params.st_num} \
#             --threshold {params.thr} \
#             --out_fig {output}
#         """


# rule fig_ST:
#     input:
#         lambda w: expand(rules.mlst.output, species=dsets[w.dset]["species"]),
#     output:
#         "results/datasets/{dset}/figs/ST_distribution.pdf",
#     conda:
#         "envs/bioinfo.yml"
#     shell:
#         """
#         python3 scripts/plot_ST.py \
#             --mlst {input} \
#             --out_fig {output}
#         """


# rule fig_metadata:
#     input:
#         rules.create_dset.output.mtd,
#     output:
#         "results/datasets/{dset}/figs/metadata.pdf",
#     conda:
#         "envs/bioinfo.yml"
#     shell:
#         """
#         python3 scripts/plot_metadata.py --metadata {input} --out_fig {output}
#         """


rule all:
    input:
        expand(rules.info_to_tsv.output, species=species),
        expand(rules.mlst.output, species=species),
        expand(rules.attotree.output.tree, species=species),
        expand(rules.assembly_to_chrom_acc.output, species=species),
