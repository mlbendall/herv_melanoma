#! /usr/bin/env python
# -*- coding: utf-8 -*-


localrules: sample_complete
rule sample_complete:
    input:
        rules.telescope.output,
        rules.kallisto.output
    output:
        touch("samples/{sampid}/completed.txt")


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


rule bowtie2_multi:
    input:
        "samples/{sampid}/unmapped.bam",
        rules.bowtie2_index.output
    output:
        "samples/{sampid}/bt2_multi.bam"
    params:
        index = config['indexes']['bowtie2']
    log:
        "samples/{sampid}/bt2_multi.log"
    conda:
        "envs/bowtie2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sampid}.XXXXXX)

picard SamToFastq\
 -I {input[0]}\
 -F $tdir/R1.fq\
 -F2 $tdir/R2.fq\
 --CLIPPING_ATTRIBUTE XT\
 --CLIPPING_ACTION 2\
 -NON_PF true\
 --TMP_DIR $tdir

(bowtie2\
 -p {threads}\
 -x {params.index}\
 -1 $tdir/R1.fq\
 -2 $tdir/R2.fq\
 -k 100\
 --very-sensitive-local\
 --score-min L,0,1.6\
 -X 1200\
 | samtools view -b > {output[0]}\
) 3>&1 1>&2 2>&3 | tee {log[0]}

chmod 600 {output[0]}
rm -rf $tdir
        '''


rule kallisto:
    input:
        "samples/{sampid}/unmapped.bam",
        rules.kallisto_index.output
    output:
        "samples/{sampid}/abundance.h5",
        "samples/{sampid}/abundance.tsv",
        "samples/{sampid}/run_info.json" 
    log:
        "samples/{sampid}/kallisto.log"
    conda:
        "envs/kallisto.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sampid}.XXXXXX)

picard SamToFastq\
 -I {input[0]}\
 -F $tdir/R1.fq\
 -F2 $tdir/R2.fq\
 --CLIPPING_ATTRIBUTE XT\
 --CLIPPING_ACTION N\
 -NON_PF true\
 --TMP_DIR $tdir

HDF5_USE_FILE_LOCKING=FALSE \
kallisto quant\
 -t {threads}\
 -b 100\
 -i {config[indexes][kallisto]}\
 -o $tdir\
 $tdir/R1.fq\
 $tdir/R2.fq\
 2>&1 | tee {log[0]}
 
mv $tdir/abundance.h5 {output[0]}
mv $tdir/abundance.tsv {output[1]}
mv $tdir/run_info.json {output[2]}

rm -rf $tdir
        '''


rule telescope:
    input:
        "samples/{sampid}/bt2_multi.bam",
        rules.telescope_annotation.output
    output:
        "samples/{sampid}/telescope.report.tsv",
        "samples/{sampid}/telescope.updated.bam"
    log:
        "samples/{sampid}/telescope.log"
    conda:
        "envs/telescope.yaml"
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sampid}.XXXXXX)

telescope assign\
 --exp_tag inform\
 --theta_prior 200000\
 --max_iter 200\
 --updated_sam\
 --outdir $tdir\
 {input[0]}\
 {input[1]}\
 2>&1 | tee {log[0]}

mv $tdir/inform-telescope_report.tsv {output[0]}
mv $tdir/inform-updated.bam {output[1]}
chmod 600 {output[1]}

rm -rf $tdir
        '''

