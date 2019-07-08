#! /usr/bin/env python
# -*- coding: utf-8 -*-

PILOT_IDS = [k for k,v in METADATA.items() if v['pilot']=='True']

localrules: symlink_abundance
rule symlink_abundance:
    input:
        "samples/{sampid}/kallisto.abundance.h5"
    output:
        temp("samples/{sampid}/abundance.h5")
    shell:
        "ln -s kallisto.abundance.h5 {output}"


rule load_sample_data:
    input:
        samples_tsv = "metadata/tcga_samples.tsv"
    output:
        "analysis/01-load_sample_data.Rdata"
    conda:
        "envs/r-tidyverse.yaml"
    script:
        "analysis/01-load_sample_data.R"

rule load_gene_data:
    input:
        ttg_tsv = config['annotations']['ttg'],
        gsym_tsv = config['annotations']['gsym'],
        herv_tsv = config['annotations']['herv_tsv'],
        l1_tsv = config['annotations']['l1_tsv']
    output:
        "analysis/01-load_gene_data.Rdata"
    conda:
        "envs/r-tidyverse.yaml"
    script:
        "analysis/01-load_gene_data.R"

rule load_tx_data:
    input:
        samp_rdata = rules.load_sample_data.output,
        gene_rdata = rules.load_gene_data.output,
        h5_files = expand("samples/{s}/abundance.h5", s=PILOT_IDS)
    output:
        "analysis/02-load_tx_data.Rdata"
    conda:
        "envs/r-dexp.yaml"
    script:
        "analysis/02-load_tx_data.R"


rule load_rtx_data:
    input:
        samp_rdata = rules.load_sample_data.output,
        gene_rdata = rules.load_gene_data.output,
        tele_files = expand("samples/{s}/telescope.report.tsv", s=PILOT_IDS)
    output:
        "analysis/03-load_rtx_data.Rdata"
    conda:
        "envs/r-dexp.yaml"
    script:
        "analysis/03-load_rtx_data.R"


rule load_metrics:
    input:
        samp_rdata = rules.load_sample_data.output,
        tele_files = expand("samples/{s}/telescope.report.tsv", s=PILOT_IDS),
        bt2_logs = expand("samples/{s}/bt2_multi.log", s=PILOT_IDS)
    output:
        "analysis/04-load_metrics.Rdata"
    conda:
        "envs/r-dexp.yaml"
    script:
        "analysis/04-load_metrics.R"

rule deseq:
    input:
        samp_rdata = rules.load_sample_data.output,
        tx_rdata = rules.load_tx_data.output,
        rtx_rdata = rules.load_rtx_data.output,
        metrics_rdata = rules.load_metrics.output
    output:
        "analysis/05-deseq.Rdata"
    conda:
        "envs/r-dexp.yaml"
    threads: snakemake.utils.available_cpu_count()
    script:
        "analysis/05-deseq.R"

rule volcano:
    input:
        deseq_rdata = rules.deseq.output,
        gene_rdata = rules.load_gene_data.output
    output:
        'analysis/06-volcano.pdf'
    conda:
        "envs/r-dexp.yaml"
    script:
        'analysis/06-volcano.R'

rule pca:
    input:
        deseq_rdata = rules.deseq.output
    output:
        'analysis/07-pca.pdf'
    conda:
        "envs/r-dexp.yaml"
    script:
        'analysis/07-pca.R'

rule gsea:
    input:
        deseq_rdata = rules.deseq.output,
        gene_rdata = rules.load_gene_data.output,
        h_gmt = 'refs/annotation/h.all.v6.2.symbols.gmt',
        c6_gmt = 'refs/annotation/c6.all.v6.2.symbols.gmt',
        c7_gmt = 'refs/annotation/c7.all.v6.2.symbols.gmt'
    output:
        pdf = 'analysis/08-gsea.pdf',
        rdata = 'analysis/08-gsea.Rdata'
    conda:
        "envs/r-dexp.yaml"
    script:
        'analysis/08-gsea.R'

rule proportions:
    input:
        samp_rdata = rules.load_sample_data.output,
        rtx_rdata = rules.load_rtx_data.output,
        metrics_rdata = rules.load_metrics.output,
        deseq_rdata = rules.deseq.output
    output:
        pdf = 'analysis/09-proportions.pdf'
    conda:
        "envs/r-dexp.yaml"
    script:
        'analysis/09-proportions.R'    

