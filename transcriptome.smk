#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule hisat2_multi:
    input:
        "samples/{sampid}/unmapped.bam",
        rules.hisat2_index.output
    output:
        "samples/{sampid}/ht2_multi.bam"
    log:
        log = "samples/{sampid}/ht2_multi.log",
        summary = "samples/{sampid}/ht2_multi.summary.log"
    conda:
        "envs/hisat2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)
        
        picard SamToFastq\
          I={input[0]} F=$tdir/R1.fq F2=$tdir/R2.fq\
          CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 NON_PF=true TMP_DIR=$tdir
        
        hisat2\
          -p {threads} --dta -k 200 --score-min L,0.0,-0.6\
          --summary-file {log.summary} --new-summary\
          --rg-id {wildcards.sampid} --rg LB:l1 --rg PL:illumina --rg PU:u1 --rg SM:s1\
          -1 $tdir/R1.fq -2 $tdir/R2.fq -x {config[indexes][hisat2]}\
        2> {log.log} | samtools view -b > {output[0]}

        chmod 600 {output[0]}
        rm -rf $tdir
        '''

####  --rna-strandness R


rule hisat2_default:
    input:
        "samples/{sampid}/unmapped.bam",
        rules.hisat2_index.output
    output:
        "samples/{sampid}/ht2_default.bam"
    log:
        log = "samples/{sampid}/ht2_default.log",
        summary = "samples/{sampid}/ht2_default.summary.log"
    conda:
        "envs/hisat2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)
        
        picard SamToFastq\
          I={input[0]} F=$tdir/R1.fq F2=$tdir/R2.fq\
          CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 NON_PF=true TMP_DIR=$tdir
        
        hisat2\
          -p {threads} --dta --score-min L,0.0,-0.6\
          --summary-file {log.summary} --new-summary\
          --rg-id {wildcards.sampid} --rg LB:l1 --rg PL:illumina --rg PU:u1 --rg SM:s1\
          -1 $tdir/R1.fq -2 $tdir/R2.fq -x {config[indexes][hisat2]}\
        2> {log.log} | samtools view -b > {output[0]}

        chmod 600 {output[0]}
        rm -rf $tdir
        '''


rule sortbam:
    input:
        "{f}.bam"
    output:
        "{f}.sorted.bam",
        "{f}.sorted.bam.bai",
    conda:
        "envs/samtools.yaml"
    threads: snakemake.utils.available_cpu_count()        
    shell:
       '''
        samtools sort -@ {threads} -T {config[tmpdir]} {input[0]} > {output[0]}
        samtools index {output[0]}
        '''


rule merge_project:
    input:
        lambda wildcards: expand("samples/{s}/ht2_multi.sorted.bam",
                                 s=BY_PROJECT['TCGA-%s' % wildcards.proj])
    output:
        'transcriptome_{proj}/merged.bam'
    conda:
        "envs/samtools.yaml"        
    threads: snakemake.utils.available_cpu_count()         
    shell:
        '''
        mkdir -p $(dirname {output})
        samtools merge -@ {threads} {output} {input}
        samtools index {output}        
        '''

rule trinity:
    input:
        'transcriptome_{proj}/merged.bam'    
    output:
        'transcriptome_{proj}/trinity_out_dir'      
    conda:
        "envs/trinity.yaml"
    threads: snakemake.utils.available_cpu_count()         
    shell:
        '''
        mkdir -p {output}
        which lfs && lfs setstripe -c 1 {output} || echo "Not lustre"
        Trinity --genome_guided_bam {input} --genome_guided_max_intron 10000\
         --max_memory 128G --CPU {threads}
        '''



'''
Trinity

--genome_guided_bam
--genome_guided_max_intron
--genome_guided_min_coverage
--genome_guided_min_reads_per_partition

Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
                --genome_guided_max_intron 10000 --CPU 6
                
Trinity --genome_guided_bam rnaseq_alignments.csorted.bam --max_memory 50G
#                --genome_guided_max_intron 10000 --CPU 6
                
 --workdir
 '''
 