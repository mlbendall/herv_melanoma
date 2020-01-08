#! /usr/bin/env python
# -*- coding: utf-8 -*-

UM_IDS = [k for k,v in METADATA.items() if v['project_id']=='TCGA-UVM']


rule load_sample_dataUM:
    input:
        samples_tsv = "metadata/tcga_samples.tsv",
    output:
        "analysisUM/01-load_sample_data.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    script:
        "analysisUM/01-load_sample_data.R"

# Sample information comes from:
# Robertson AG, Shih J, Yau C, et al. Integrative Analysis Identifies Four 
# Molecular and Clinical Subsets in Uveal Melanoma. Cancer Cell. 
# 2017;32(2):204-220.e15. doi:10.1016/j.ccell.2017.07.003
#
# Table S1 is saved as "metadata/mmc2.xlsx"

rule load_clin_dataUM:
    input:
        samp_rdata = rules.load_sample_dataUM.output,
        tableS1 = "metadata/mmc2.xlsx"
    output:
        "analysisUM/01-load_clin_data.Rdata"    
    conda:
        "envs/r-analysisUM.yaml"
    script:
        "analysisUM/01-load_clin_data.R"


rule load_tx_dataUM:
    input:
        samp_rdata = rules.load_sample_dataUM.output,
        gene_rdata = rules.load_gene_data.output,
        h5_files = expand("samples/{s}/abundance.h5", s=UM_IDS)
    output:
        "analysisUM/02-load_tx_data.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    script:
        "analysis/02-load_tx_data.R"


rule load_rtx_dataUM:
    input:
        samp_rdata = rules.load_sample_dataUM.output,
        gene_rdata = rules.load_gene_data.output,
        tele_files = expand("samples/{s}/telescope.report.tsv", s=UM_IDS)
    output:
        "analysisUM/03-load_rtx_data.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    script:
        "analysis/03-load_rtx_data.R"


rule load_metricsUM:
    input:
        samp_rdata = rules.load_sample_dataUM.output,
        tele_files = expand("samples/{s}/telescope.report.tsv", s=UM_IDS),
        bt2_logs = expand("samples/{s}/bt2_multi.log", s=UM_IDS)
    output:
        "analysisUM/04-load_metrics.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    script:
        "analysis/04-load_metrics.R"

rule deseq_intUM:
    input:
        samp_rdata = rules.load_sample_dataUM.output,
        tx_rdata = rules.load_tx_dataUM.output,
        rtx_rdata = rules.load_rtx_dataUM.output,
        metrics_rdata = rules.load_metricsUM.output
    output:
        "analysisUM/05-deseq_int.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    threads: snakemake.utils.available_cpu_count()
    script:
        "analysisUM/05-deseq_int.R"


rule deseq_SCNA:
    input:
        clin_rdata = rules.load_clin_dataUM.output,
        tx_rdata = rules.load_tx_dataUM.output,
        rtx_rdata = rules.load_rtx_dataUM.output,
        metrics_rdata = rules.load_metricsUM.output
    output:
        "analysisUM/05-deseq_SCNA.Rdata"
    conda:
        "envs/r-analysisUM.yaml"
    threads: snakemake.utils.available_cpu_count()
    script:
        "analysisUM/05-deseq_SCNA.R"

 
rule pcaUM:
    input:
        clin_rdata = rules.load_clin_dataUM.output,    
        deseq_rdata = rules.deseq_intUM.output
    output:
        'analysisUM/07-pca.pdf'
    conda:
        "envs/r-analysisUM.yaml"
    script:
        'analysisUM/07-pca.R'

"""
rule um_analysis:
    input:
        rules.deseq_intUM.output,
"""