# Melanoma Retrotranscriptome

## Sample metadata

```bash
scripts/gdcquery.py TCGA-UVM > metadata/TCGA-UVM.tsv
scripts/gdcquery.py TCGA-SKCM > metadata/TCGA-SKCM.tsv
```

```python
scripts/make_sample_table.py
```

## References

```bash
snakemake bowtie2_index -n --profile ./wcm.profile
snakemake kallisto_index -n --profile ./wcm.profile
```




```bash
mkdir -p refs
mkdir -p snakelogs
```

```bash
snakemake -j 20 --use-conda -pr all_references
```


```bash
mkdir -p refs
snakemake -j 20 --use-conda -pr all_references
```


```bash
snakemake --profile ./c1.profile samples/TCGA-DA-A95X-06A/completed.txt
```

```bash
grep 'TCGA-UVM' metadata/tcga_samples.tsv
mkdir -p transcriptome_uvm
```


```bash
snakemake -j 20 --use-conda transcriptome_UVM/merged.bam
```

`pybedtools.bedtool.BedTool.genome_coverage(…)`

```bash
bedtools genomecov -ibam samples/TCGA-YZ-A984-01A/ht2_multi.sorted.bam -bga > tmp.bed
bedtools genomecov -ibam samples/TCGA-YZ-A985-01A/ht2_multi.sorted.bam -bga > tmp2.bed
bedtools genomecov -ibam samples/TCGA-YZ-A983-01A/ht2_multi.sorted.bam -bga > tmp3.bed
bedtools unionbedg -empty -g genome.file -i tmp.bed tmp2.bed tmp3.bed > union.txt


cat tmp.bed | awk '$4 > 9' | bedtools merge -d 1000 -c 4 -o median > tmp1.m.bg
cat tmp2.bed | awk '$4 > 9' | bedtools merge -d 1000 -c 4 -o median > tmp2.m.bg
cat tmp3.bed | awk '$4 > 9' | bedtools merge -d 1000 -c 4 -o median > tmp3.m.bg
bedtools unionbedg -empty -g genome.file -i tmp1.m.bg tmp2.m.bg tmp3.m.bg > union.txt

python -c 'import sys; z = (l.strip().split('\t') for l in sys.stdin); print('\n'.join('\t'.join([l[0],l[1],l[2]]) for l in z))'


#!/usr/bin/env python
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import sys

for l in sys.stdin:
    f = l.strip().split('\t')
    total_cov = sum(map(float, f[3:]))
    if total_cov <= 0:
        print('%s\t%s\t%s\t%g' % (f[0], f[1], f[2], total_cov), file=sys.stdout)



```