set -euo pipefail

# Complie

gcc -O2 -Wall -Wextra -o locusDefiner locusDefiner.c -lz -lm

#  Downlaod Example File

if [ ! -e ./example/QTL_BMI.EUR.gwama.sumstats.txt.gz ]; then

  wget -O ./example/QTL_BMI.EUR.gwama.sumstats.txt.gz \
    https://g-fce312.fd635.8443.data.globus.org/sumstats_downsized/EUR/QTL_BMI.EUR.gwama.sumstats.txt.gz

fi

# Execution

time ./locusDefiner \
  --input ./example/QTL_BMI.EUR.gwama.sumstats.txt.gz \
  --output ./example_loci.tsv \
  --thr 5e-8 \
  --window 500000

n=$(cat ./example_loci.tsv | sed '1d' | wc -l)
echo ${n} loci detected # 707 loci
