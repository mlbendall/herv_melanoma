#! /usr/bin/env python

from collections import Counter
import random

rowdict = []

# Parse SKCM
tsv = (l.strip('\n').split('\t') for l in open('metadata/TCGA-SKCM.tsv'))
header = next(tsv)
for row in tsv:
    rowdict.append(dict(zip(header,row)))

# Parse UVM
tsv = (l.strip('\n').split('\t') for l in open('metadata/TCGA-UVM.tsv'))
header = next(tsv)
for row in tsv:
    rowdict.append(dict(zip(header,row)))

# Identify cases with multiple samples
c = Counter(d['cases.0.submitter_id'] for d in rowdict)
multiple = set([k for k,v in c.most_common() if v>1])

# Construct new table
newtable = []
for rowd in rowdict:
    assert rowd['id'] == rowd['file_id']
    newrow = [
        rowd['cases.0.project.project_id'],
        rowd['cases.0.submitter_id'],
        rowd['cases.0.samples.0.submitter_id'],
        rowd['cases.0.samples.0.sample_type'],
        rowd['id'],
        rowd['file_name'],
        rowd['cases.0.case_id'],
        rowd['md5sum'],
        'False',
    ]
    newtable.append(newrow)

random.seed(12345)
selected = []
selected += random.sample([i for i,l in enumerate(newtable) if l[0] == 'TCGA-UVM' and l[1] not in multiple], 50)
selected += random.sample([i for i,l in enumerate(newtable) if l[0] == 'TCGA-SKCM' and l[1] not in multiple and l[3]=='Primary Tumor'], 50)
selected += random.sample([i for i,l in enumerate(newtable) if l[0] == 'TCGA-SKCM' and l[1] not in multiple and l[3]=='Metastatic'], 50)

assert not any(newtable[idx][1] in multiple for idx in selected), 'Samples selected have multiple'

for idx in selected:
    newtable[idx][-1] = 'True'

newheader = ['project_id', 'case_id', 'sample_id', 'sample_type', 'file_uuid', 'file_name', 'case_uuid', 'md5sum', 'pilot']
with open('metadata/tcga_samples.tsv', 'w') as outh:
    print('\t'.join(newheader), file=outh)
    print('\n'.join('\t'.join(r) for r in newtable), file=outh)
