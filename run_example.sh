set -euo pipefail

gcc -O2 -Wall -Wextra -o locusDefiner locusDefiner.c -lz -lm

./locusDefiner \
  --input ./example/QTL_BMI.EUR.gwama.sumstats.txt.gz \
  --output ./example_loci.tsv \
  --thr 5e-8 \
  --window 500000

n=$(cat ./example_loci.tsv | sed '1d' | wc -l)
echo ${n} loci detected
