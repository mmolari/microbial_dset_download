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


rule mash_triangle:
    input:
        rules.chromosome_fa.output,
    output:
        "results/{species}/mash_triangle.tsv",
    params:
        k=21,  # k-mer size
        s=50000,  # sketch size
        cores=4,  # number of threads
    conda:
        "envs/mash.yml"
    shell:
        """
        mash triangle \
            -k {params.k} \
            -s {params.s} \
            -p {params.cores} \
            {input}/*.fa > {output}
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
        s=50000,  # sketch size
        fof="data/species/{species}/attotree.fof",  # temporary file of filenames
    shell:
        """
        # 1. Build the file-of-filenames once.
        find {input.fas} -maxdepth 1 -name '*.fa' -print > {params.fof}

        # 2. Run attotree in list-input mode.
        attotree -L {params.fof} \
                 -o {output.tree} \
                 -t {threads} \
                 -k {params.k} \
                 -s {params.s}

        # 3. Clean up the list file.
        rm -f {params.fof}
        """


rule refine_tree:
    input:
        tree="data/species/{species}/attotree.nwk",
    output:
        "results/{species}/attotree.nwk",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/refine_tree.py \
            --tree {input.tree} \
            --out_tree {output}
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


rule combine_metadata:
    input:
        metadata="data/species/{species}/info.tsv",
        chrom="data/species/{species}/assembly_to_chrom.tsv",
        mlst="data/species/{species}/mlst.tsv",
    output:
        "results/{species}/combined_metadata.csv",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/combine_metadata.py \
            --metadata_file {input.metadata} \
            --chromosome_mapping_file {input.chrom} \
            --mlst_file {input.mlst} \
            --output {output}
        """


rule plot_metadata_overview:
    input:
        rules.combine_metadata.output,
    output:
        "results/{species}/metadata_overview.png",
    params:
        species_name=lambda w: species[w.species],
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/plot_metadata_overview.py \
            --df {input} \
            --species '{params.species_name}' \
            --fig {output}
        """


# Plot tree with metadata using plot_tree.py
rule plot_tree:
    input:
        tree=rules.refine_tree.output,
        metadata=rules.combine_metadata.output,
    output:
        "results/{species}/tree_metadata.png",
    conda:
        "envs/bioinfo.yml"
    shell:
        """
        python3 scripts/plot_tree.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --fig {output}
        """


rule clean:
    shell:
        """
        rm -f results/*/metadata_overview.png
        """


rule all:
    input:
        expand(rules.refine_tree.output, species=species),
        expand(rules.combine_metadata.output, species=species),
        expand(rules.plot_metadata_overview.output, species=species),
        expand(rules.plot_tree.output, species=species),
        expand(rules.mash_triangle.output, species=species),
