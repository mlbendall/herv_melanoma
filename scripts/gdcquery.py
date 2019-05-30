#! /usr/bin/env python

import requests
import json
from collections import defaultdict
import sys

def query_rnaseq(project):
    endpoint = 'https://api.gdc.cancer.gov/files'
    fields = [
        'experimental_strategy',
        'cases.project.project_id',
        'file_id',
        'file_name',
        'md5sum',
        'cases.submitter_id',
        'cases.case_id',
        'data_category',
        'data_type',
        'cases.samples.tumor_descriptor',
        'cases.samples.tissue_type',
        'cases.samples.sample_type',
        'cases.samples.sample_type_id',
        'cases.samples.submitter_id',
        'cases.samples.sample_id',
        'cases.samples.portions.analytes.aliquots.aliquot_id',
        'cases.samples.portions.analytes.aliquots.submitter_id',
    ]
    filters = {
        'op': 'and',
        'content': [
            {"op":"=", "content": {"field":"experimental_strategy", "value":"RNA-Seq"}},
            {"op":"=", "content": {"field":"data_category", "value":"Sequencing Reads"}},
            {"op":"=", "content": {"field":"data_format", "value":"BAM"}},
            {"op":"=", "content": {"field":"cases.project.project_id", "value": project}},
        ]
    }
    query = {
        'format': 'TSV',
        'fields': ','.join(fields),
        'size':'2000',
        'filters': filters
    }
    headers = {"Content-Type": "application/json"}
    response = requests.post(endpoint, headers=headers, json=query)
    return [r.split('\t') for r in response.text.strip('\r\n').split('\r\n')]

if __name__ == '__main__':
    if not len(sys.argv) == 2:
       print('USAGE: %s <project_id>' % sys.argv[0], file=sys.stderr)
       sys.exit(1)
    tsv = query_rnaseq(sys.argv[1])
    print('\n'.join('\t'.join(row) for row in tsv), file=sys.stdout)
