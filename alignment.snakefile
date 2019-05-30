#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule bowtie2_multi:
    input:
        "samples/{sampid}/unmapped.bam",
        ancient(expand(config['indexes']['bowtie2'] + ".{i}.bt2", i = range(1, 5))),
        ancient(expand(config['indexes']['bowtie2'] + ".rev.{i}.bt2", i = range(1, 3)))        
    output:
        "samples/{sampid}/bt2_multi.bam"
    log:
        "samples/{sampid}/bt2_multi.log"
    conda:
        "envs/bowtie2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        'tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)'
        ' && '
        'picard'
        ' SamToFastq'
        ' I={input[0]}'
        ' F=$tdir/R1.fq F2=$tdir/R2.fq'
        ' CLIPPING_ATTRIBUTE=XT'
        ' CLIPPING_ACTION=2'
        ' NON_PF=true'
        ' TMP_DIR=$tdir'
        ' && '
        '(bowtie2'
        ' -p {threads}'
        ' -x {config[indexes][bowtie2]}'
        ' -1 $tdir/R1.fq -2 $tdir/R2.fq'
        ' -k 100'
        ' --very-sensitive-local'
        ' --score-min L,0,1.6'
        ' -X 1200'
        ' | '
        'samtools view -b > {output[0]}'
        ') 3>&1 1>&2 2>&3 | tee {log[0]}'
        ' && '
        'chmod 600 {output[0]}'
        ' && '
        'rm -rf $tdir'


rule kallisto:
    input:
        "samples/{sampid}/unmapped.bam"
    output:
        "samples/{sampid}/kallisto.abundance.h5",
        "samples/{sampid}/kallisto.abundance.tsv",
        "samples/{sampid}/kallisto.run_info.json" 
    log:
        "samples/{sampid}/kallisto.log"
    conda:
        "envs/kallisto.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        'tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)'
        ' && '
        'picard'
        ' SamToFastq'
        ' I={input}'
        ' F=$tdir/R1.fq F2=$tdir/R2.fq'
        ' CLIPPING_ATTRIBUTE=XT'
        ' CLIPPING_ACTION=N'
        ' NON_PF=true'
        ' TMP_DIR=$tdir'
        ' && '        
        'HDF5_USE_FILE_LOCKING=FALSE '
        'kallisto'
        ' quant'
        ' -t {threads}'
        ' -b 100'
        ' -i {config[indexes][kallisto]}'
        ' -o $tdir'
        ' $tdir/R1.fq $tdir/R2.fq'
        ' 2>&1 | tee {log[0]}'
        ' && '
        'mv $tdir/abundance.h5 {output[0]} && '
        'mv $tdir/abundance.tsv {output[1]} && '
        'mv $tdir/run_info.json {output[2]}'
        ' && '
        'rm -rf $tdir'


rule telescope:
    input:
        aln = "samples/{sampid}/bt2_multi.bam",
        ann = config['annotations']['retro']
    output:
        "samples/{sampid}/telescope.report.tsv",
        "samples/{sampid}/telescope.updated.bam"
    log:
        "samples/{sampid}/telescope.log"
    conda:
        "envs/telescope.yaml"
    shell:
        'tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)'
        ' && '
        'telescope'
        ' assign'
        ' --exp_tag inform'
        ' --theta_prior 200000'
        ' --max_iter 200'
        ' --updated_sam'
        ' --outdir $tdir'
        ' {input[0]}'
        ' {input[1]}'
        ' 2>&1 | tee {log[0]}'
        ' && '
        'mv $tdir/inform-telescope_report.tsv {output[0]} &&'
        'mv $tdir/inform-updated.bam {output[1]} &&'
        'chmod 600 {output[1]}'
        ' && '
        'rm -rf $tdir'
