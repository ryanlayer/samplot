This directory holds test data for samplot. The files have been sampled to contain only reads for the following regions:
```
4       115922418       115938188       DEL
X       101031678       101090808       DUP
```

These regions contain the following variants:
```
4       115928726       115931880       DEL
X       101055330       101067156       DUP
```

The VCF file `NA12878.trio.svt.subset.vcf` has the above SVs and can be used for testing `samplot_vcf.sh` or SV-plaudit's `annotate.py`.
