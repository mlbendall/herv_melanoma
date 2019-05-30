#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""
Download data from GDC
"""
localrules: gdc_download, gdc_download_check

rule gdc_download:
    output:
        "samples/{sampid}/original.bam"
    params:
        uuid = lambda wildcards: METADATA[wildcards.sampid]['file_uuid']
    shell:
        'mkdir -p $(dirname {output[0]})'
        ' && '
        'curl'
        ' --header "X-Auth-Token: $(<{config[gdc_token_file]})"'
        ' https://api.gdc.cancer.gov/data/{params.uuid}'
        ' > {output[0]}'
        ' && '
        'chmod 600 {output[0]}'

rule gdc_download_check:
    input:
        "samples/{sampid}/original.bam"
    output:
        touch("samples/{sampid}/download.ok")
    params:
        md5sum = lambda wildcards: METADATA[wildcards.sampid]['md5sum']
    shell:
        'echo {params.md5sum} {input} | md5sum -c -'

"""
Revert mapped data to unmapped
"""
# Default SAM attributes cleared by RevertSam
attr_revertsam = ['NM', 'UQ', 'PG', 'MD', 'MQ', 'SA', 'MC', 'AS']
# SAM attributes output by STAR
attr_star = ['NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS', 'uT']
# Additional attributes to clear
ALN_ATTRIBUTES = list(set(attr_star) - set(attr_revertsam))

rule revert_and_mark_adapters:
    input:
        "samples/{sampid}/original.bam",
        "samples/{sampid}/download.ok"
    output:
        "samples/{sampid}/unmapped.bam"
    log:
        "samples/{sampid}/revert_bam.log",
        "samples/{sampid}/mark_adapters.log",
        "samples/{sampid}/mark_adapters.metrics"
    conda:
        "envs/picard.yaml"
    params:
        attr_to_clear = expand("ATTRIBUTE_TO_CLEAR={a}", a=ALN_ATTRIBUTES)
    shell:
        'picard'
        ' RevertSam'
        ' I={input[0]}'
        ' O=/dev/stdout'
        ' SANITIZE=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT'
        ' {params.attr_to_clear}'
        ' TMP_DIR={config[tmpdir]}'
        ' 2> {log[0]}'
        ' | '
        'picard'
        ' MarkIlluminaAdapters'
        ' I=/dev/stdin'
        ' O={output[0]}'
        ' M={log[2]}'
        ' COMPRESSION_LEVEL=5'
        ' TMP_DIR={config[tmpdir]}'
        ' 2> {log[1]}'
        ' && '
        ' chmod 600 {output[0]}'


# rule revert_bam:
#     input:
#         "samples/{sampid}/original.bam",
#         "samples/{sampid}/download.ok"
#     output:
#         temp("samples/{sampid}/unmapped.tmp.bam")
#     log:
#         "samples/{sampid}/revert_bam.log"
#     conda:
#         "envs/picard.yaml" 
#     params:
#         attr_to_clear = expand("ATTRIBUTE_TO_CLEAR={a}", a=ALN_ATTRIBUTES)
#     shell:
#         'picard'
#         ' RevertSam'
#         ' I={input[0]}'
#         ' O={output[0]}'
#         ' SANITIZE=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT'
#         ' {params.attr_to_clear}'
#         ' TMP_DIR={config[tmpdir]}'
#         ' 2> {log[0]}'
# 
# 
# rule mark_adapters:
#     input:
#         rules.revert_bam.output
#     output:
#         "samples/{sampid}/unmapped.bam"
#     log:
#         "samples/{sampid}/mark_adapters.log",
#         "samples/{sampid}/mark_adapters.metrics"
#     conda:
#         "envs/picard.yaml"
#     shell:
#         'picard'
#         ' MarkIlluminaAdapters'
#         ' I={input[0]}'
#         ' O={output[0]}'
#         ' M={log[1]}'
#         ' COMPRESSION_LEVEL=5'
#         ' TMP_DIR={config[tmpdir]}'
#         ' 2> {log[0]}'
#         ' && '
#         ' chmod 600 {output[0]}'
# 
# 
# rule revert_and_mark_adapters:
#     input:
#         rules.revert_bam.input
#     output:
#         rules.mark_adapters.output,
#         touch("samples/{sampid}/rma.ok")
# 
# rule ubam_to_fastq:
#     input:
#         "samples/{sampid}/unmapped.bam"
#     output:
#     shell:
#         'picard'
#         ' SamToFastq'
#         ' I={input}'
#         ' F='
#         ' F2'
#         ' CLIPPING_ATTRIBUTE=XT'
#         ' CLIPPING_ACTION=2'
#         ' TMP_DIR={config[tmpdir]}'
# 
# 
# 


