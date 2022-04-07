#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule gdc_download:
    """ Download BAM from Genomic Data Commons (GDC)
    """
    input:
        config['gdc_token_file']
    output:
        temp("samples/{sampid}/original.bam")
    params:
        uuid = lambda wc: METADATA[wc.sampid]['file_uuid'],
        md5sum = lambda wc: METADATA[wc.sampid]['md5sum']
    shell:
        '''
mkdir -p $(dirname {output[0]})

curl\
 -H "X-Auth-Token: $(<{input[0]})"\
 https://api.gdc.cancer.gov/data/{params.uuid}\
 > {output[0]}

echo {params.md5sum} {output[0]} | md5sum -c -
chmod 600 {output[0]}
        '''


# Default SAM attributes cleared by RevertSam
attr_revertsam = ['NM', 'UQ', 'PG', 'MD', 'MQ', 'SA', 'MC', 'AS']
# SAM attributes output by STAR
attr_star = ['NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS', 'uT']
# Additional attributes to clear
ALN_ATTRIBUTES = list(set(attr_star) - set(attr_revertsam))

rule revert_and_mark_adapters:
    """ Create unmapped BAM (uBAM) from aligned BAM
    """
    input:
        "samples/{sampid}/original.bam"
    output:
        protected("samples/{sampid}/unmapped.bam")
    log:
        "samples/{sampid}/revert_bam.log",
        "samples/{sampid}/mark_adapters.log",
        "samples/{sampid}/mark_adapters.metrics"
    conda:
        "envs/utils.yaml"
    params:
        attr_to_clear = expand("--ATTRIBUTE_TO_CLEAR {a}", a=ALN_ATTRIBUTES),
        tmpdir = config['local_tmp']
    shell:
        '''
picard RevertSam\
 -I {input[0]}\
 -O /dev/stdout\
 --SANITIZE true\
 --COMPRESSION_LEVEL 0\
 --VALIDATION_STRINGENCY SILENT\
 {params.attr_to_clear}\
 --TMP_DIR {params.tmpdir}\
 2> {log[0]}\
 | picard MarkIlluminaAdapters\
 -I /dev/stdin\
 -O {output[0]}\
 -M {log[2]}\
 --COMPRESSION_LEVEL 5\
 --TMP_DIR {params.tmpdir}\
 2> {log[1]}

chmod 600 {output[0]}
        '''
