#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()

from collections import defaultdict

# TOKEN = 'gdc-user-token.2019-05-13T21_17_16.257Z.txt'

wildcard_constraints:
    caseid = 'TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}',
    sampid = 'TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-\d{2}[A-Z]',
    uuid = '[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-4[0-9A-Fa-f]{3}-[89ABab][0-9A-Fa-f]{3}-[0-9A-Fa-f]{12}'

configfile: "config.yaml"

# Parse metadata
METADATA = {}
tsv = (l.strip('\n').split('\t') for l in open('metadata/tcga_samples.tsv', 'r'))
header = next(tsv)
for row in tsv:
    d = dict(zip(header, row))
    METADATA[d['sample_id']] = d

BY_PROJECT = defaultdict(list)
for k,v in METADATA.items():
    BY_PROJECT[v['project_id']].append(k)

localrules: pilot, all, sample_complete

rule all:
    input:
        expand("samples/{s}/completed.txt", s=METADATA.keys())

rule pilot:
    input:
        expand("samples/{s}/completed.txt", s=[k for k,v in METADATA.items() if v['pilot']=='True'])
        #expand("samples/{s}/completed.txt", s=['TCGA-DA-A95X-06A', 'TCGA-EB-A57M-01A', 'TCGA-EB-A6R0-01A', ])
    output:
        "pilot_results.tgz"        
    shell:
        'mkdir -p pilot_results && '
        'fstr="{input}"'
        ' && '
        'for f in $fstr; do'
        ' s=$(basename $(dirname $f)) &&'
        ' mkdir -p pilot_results/$s &&'
        ' cp samples/$s/bt2_multi.log pilot_results/$s &&'
        ' cp samples/$s/kallisto.abundance.h5 pilot_results/$s &&'
        ' cp samples/$s/kallisto.abundance.tsv pilot_results/$s &&'
        ' cp samples/$s/kallisto.run_info.json pilot_results/$s &&'
        ' cp samples/$s/telescope.report.tsv pilot_results/$s; '
        'done'
        ' && '
        'tar -czf {output} pilot_results'

''' References '''
include: "refs.snakefile"

localrules: all_references

rule all_references:
    input:
        rules.kallisto_index.output,
        rules.bowtie2_index.output,
        rules.telescope_annotation.output

''' '''
include: "utils.snakefile"

''' '''
include: "alignment.snakefile"


localrules: sample_complete

rule sample_complete:
    input:
        rules.kallisto.output,
        rules.telescope.output
    output:
        touch("samples/{sampid}/completed.txt")
    shell:
        'rm -f samples/{wildcards.sampid}/original.bam'


include: "analysis.snakefile"
include: "analysisUM.snakefile"

include: "transcriptome.smk"

