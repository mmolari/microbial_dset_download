configfile: "config.yml"

ref_dict = config["ref_acc"]

rule download_ref:
    output:
        "data/ref/{r_acc}.fa"
    conda:
        "envs/ncbi.yml"
    shell:
        """
        datasets download genome accession {wildcards.r_acc} \
            --filename {wildcards.r_acc}.zip
        unzip {wildcards.r_acc}.zip -d {wildcards.r_acc}
        mv {wildcards.r_acc}/ncbi_dataset/data/*/*.fna {output}
        rm -r {wildcards.r_acc} {wildcards.r_acc}.zip
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
        temp("data/dset.zip")
    conda:
        "envs/ncbi.yml"
    shell:
        """
        datasets download genome taxon 'staphylococcus aureus' \
            --assembly-level complete \
            --assembly-source 'RefSeq' \
            --exclude-atypical \
            --exclude-multi-isolate \
            --filename {output}
        """

rule expand_dataset:
    input:
        rules.download_dataset.output
    output:
        ds=directory("data/dset"),
        info=temp("data/dset_info.jsonl")
    conda:
        "envs/ncbi.yml"
    shell:
        """
        unzip {input} -d tmp
        mkdir -p {output.ds}
        for f in tmp/ncbi_dataset/data/*/*.fna; do
            mv $f {output.ds}/$(basename $f)
        done

        mv tmp/ncbi_dataset/data/*.jsonl {output.info}
        rm -r tmp
        """

rule info_to_tsv:
    input:
        rules.expand_dataset.output.info
    output:
        "data/dset_info.tsv"
    conda:
        "envs/ncbi.yml"
    shell:
        """
        dataformat tsv genome --inputfile {input} > {output}
        """

rule all:
    input:
        # rules.download_summary.output,
        rules.expand_dataset.output,
        rules.info_to_tsv.output,
        expand("data/ref/{r_acc}.fa", r_acc=ref_dict.values())
