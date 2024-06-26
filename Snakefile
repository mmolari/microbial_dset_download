configfile: "config.yml"


dsets = config["dsets"]
species = config["species"]


rule download_ref:
    output:
        "data/ref/record/{ref_acc}.fa",
    conda:
        "envs/ncbi.yml"
    shell:
        """
        datasets download genome accession {wildcards.ref_acc} \
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


# rule download_summary:
#     output:
#         "data/summary_info.json"
#     conda:
#         "envs/ncbi.yml"
#     shell:
#         """
#         datasets summary genome taxon 'staphylococcus aureus' \
#             --assembly-level complete \
#             --assembly-source 'RefSeq' \
#             --exclude-atypical \
#             --exclude-multi-isolate \
#             > {output}
#         """


rule download_dataset:
    output:
        "data/species/{species}/ncbi.zip",
    params:
        ncbi_sp=lambda w: species[w.species],
    conda:
        "envs/ncbi.yml"
    shell:
        """
        datasets download genome taxon '{params.ncbi_sp}' \
            --assembly-level complete \
            --assembly-source 'RefSeq' \
            --exclude-atypical \
            --exclude-multi-isolate \
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
    conda:
        "envs/ncbi.yml"
    shell:
        """
        dataformat tsv genome --inputfile {input} > {output}
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


rule mash_dist:
    input:
        fas=rules.chromosome_fa.output,
        ref=rules.ref_to_chromosome.output,
    output:
        "data/species/{species}/mash_dist/{ref_acc}.tsv",
    conda:
        "envs/mash.yml"
    shell:
        """
        mash dist -s 10000 {input.ref} {input.fas}/*.fa > {output}
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
        mlst --scheme {params.scheme} {input}/*.fa > {output}
        """


rule fig_mash_dist:
    input:
        mash=lambda w: expand(
            rules.mash_dist.output,
            species=dsets[w.dset]["species"],
            ref_acc=dsets[w.dset]["ref_acc"],
        ),
        mlst=lambda w: expand(
            rules.mlst.output,
            species=dsets[w.dset]["species"],
        ),
    output:
        "results/{dset}/figs/mash_dist.pdf",
    params:
        st_num=lambda w: dsets[w.dset]["strain"],
        thr=lambda w: dsets[w.dset]["threshold"],
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/plot_mash_dist.py \
            --mash {input.mash} \
            --mlst {input.mlst} \
            --ST_num {params.st_num} \
            --threshold {params.thr} \
            --out_fig {output}
        """


rule fig_ST:
    input:
        lambda w: expand(rules.mlst.output, species=dsets[w.dset]["species"]),
    output:
        "results/{dset}/figs/ST_distribution.pdf",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/plot_ST.py \
            --mlst {input} \
            --out_fig {output}
        """


rule create_dset:
    input:
        mlst=lambda w: expand(rules.mlst.output, species=dsets[w.dset]["species"]),
        mash=lambda w: expand(
            rules.mash_dist.output,
            species=dsets[w.dset]["species"],
            ref_acc=dsets[w.dset]["ref_acc"],
        ),
        mtd=lambda w: expand(rules.info_to_tsv.output, species=dsets[w.dset]["species"]),
        acm=lambda w: expand(
            rules.assembly_to_chrom_acc.output, species=dsets[w.dset]["species"]
        ),
        fas=lambda w: expand(
            rules.chromosome_fa.output, species=dsets[w.dset]["species"]
        ),
    output:
        seqs=directory("results/datasets/{dset}/fa"),
        mtd="results/datasets/{dset}/metadata.tsv",
        mlst="results/datasets/{dset}/mlst.tsv",
    params:
        thr=lambda w: dsets[w.dset]["threshold"],
        st=lambda w: dsets[w.dset]["strain"],
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/create_dset.py \
            --mlst {input.mlst} \
            --mash {input.mash} \
            --mtd {input.mtd} \
            --acc_map {input.acm} \
            --fas {input.fas} \
            --ST {params.st} \
            --threshold {params.thr} \
            --out_dir {output.seqs} \
            --out_mtd {output.mtd} \
            --out_mlst {output.mlst}
        """


rule fig_metdata:
    input:
        rules.create_dset.output.mtd,
    output:
        "results/{dset}/figs/metadata.pdf",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/plot_metadata.py --metadata {input} --out_fig {output}
        """


rule all:
    input:
        expand(rules.fig_ST.output, dset=dsets.keys()),
        expand(rules.fig_mash_dist.output, dset=dsets.keys()),
        expand(rules.create_dset.output, dset=dsets.keys()),
        expand(rules.fig_metdata.output, dset=dsets.keys()),
