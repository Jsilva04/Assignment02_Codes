set -e
sed '1d' metadata.csv | while read meta
do
  accn=$(echo $meta | cut -d "," -f 1)
  name=$(echo $meta | cut -d "," -f 39)
  fasterq-dump ${accn}
  head -n 200000 ${accn}_1.fastq | gzip > ${name}_1.fastq.gz
  rm ${accn}_1.fastq
  head -n 200000 ${accn}_2.fastq | gzip > ${name}_2.fastq.gz
  rm ${accn}_2.fastq
done
