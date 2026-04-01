# Locus Definer
Define a locus from GWAS summary statistics
``` exemple
locusDefiner \
  --input YOUR_GWAS.sumstats.txt.gz \
  --output YOUR_LOCI.tsv \
  --thr 5e-8 \
  --window 500000
```

options
```
--col-chr NAME [default CHR]
--col-pos NAME [default POS]
--col-p NAME [default P]
--p-scale (p|neglog10) [default p]
```
