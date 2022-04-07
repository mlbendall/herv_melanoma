#! /usr/bin/env python
# -*- coding: utf-8 -*-

UM_IDS = [k for k,v in METADATA.items() if v['project_id']=='TCGA-UVM']


rule load_sample_data:
    input:
        samples_tsv = "metadata/tcga_samples.tsv",
    output:
        "analysis2/01-sample_data.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/01-sample_data.R"


rule load_clinical_data:
    input:
        samp_rdata = rules.load_sample_data.output,
        tableS1 = "metadata/mmc2.xlsx",
        cibersort_table = "metadata/TCGA.Kallisto.fullIDs.cibersort.relative.tsv"
    output:
        "analysis2/02-clinical_data.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/02-clinical_data.R"


rule build_ens_DB:
    input:
        "refs/downloads/Homo_sapiens.GRCh38.99.gtf.gz"
    output:
        'refs/annotation/Homo_sapiens.GRCh38.99.sqlite'
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/03a-build_ens_DB.R"


rule load_gene_data:
    input:
        ttg_tsv = config['annotations']['ttg'],
        gsym_tsv = config['annotations']['gsym'],
        herv_tsv = config['annotations']['herv_tsv'],
        l1_tsv = config['annotations']['l1_tsv'],
        tedb_tsv = 'refs/annotation/TE_annotation.v2.0.tsv',
        geneinfo_tsv = config['annotations']['gene_info']
    output:
        "analysis2/03-gene_data.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/03-gene_data.R"


rule load_count_data:
    input:
        samp_rdata = rules.load_sample_data.output,
        gene_rdata = rules.load_gene_data.output,
        h5_files = expand("samples/{s}/abundance.h5", s=UM_IDS),
        tele_files = expand("samples/{s}/telescope.report.tsv", s=UM_IDS)        
    output:
        "analysis2/04-count_data.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/04-count_data.R"


rule filter_count_data:
    input:
        count_rdata = rules.load_count_data.output
    output:
        "analysis2/05-filter_counts.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/05-filter_counts.R"


rule unsupervised:
    input:
        samp_rdata = rules.load_sample_data.output,
        clin_rdata = rules.load_clinical_data.output,
        filt_rdata = rules.filter_count_data.output
    output:
        rdata = "analysis2/06-unsupervised.Rdata",
        outdir = directory("analysis2/06-unsupervised")
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/06-unsupervised.R"


rule dendro:
    input:
        samp_rdata = rules.load_sample_data.output,
        clin_rdata = rules.load_clinical_data.output,
        clust_rdata = rules.unsupervised.output.rdata
    output:
        outpdf = "analysis2/06-unsupervised/colorbars.pdf",
        rdata = "analysis2/06a-dendro.Rdata"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/06a-dendro.R"


rule survival:
    input:
        samp_rdata = rules.load_sample_data.output,
        clin_rdata = rules.load_clinical_data.output,
        clust_rdata = rules.unsupervised.output.rdata
    output:
        outpdf = "analysis2/06-unsupervised/survival.pdf"
    conda:
        "envs/r-analysis2.yaml"
    script:
        "analysis2/06b-survival.R"



