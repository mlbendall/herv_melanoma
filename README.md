```bash
mkdir -p refs
mkdir -p snakelogs

```
snakemake -j 20 --use-conda -pr all_references
```

```bash
mkdir -p refs
snakemake -j 20 --use-conda -pr all_references
```

```bash
snakemake --profile ./c1.profile samples/TCGA-DA-A95X-06A/completed.txt
```