DIR="/path/to/cellranger_output/"
ID=$1
OUTPUT=2
cd ${OUTPUT}
mkdir -p fq bam vcf cts
python renamer.py --bam ${DIR}/${ID}/outs/gex_possorted_bam.bam --barcodes ${DIR}/${ID}/outs/barcodes.tsv --out fq/${ID}_rename.fq
minimap2 -ax splice -t 8 -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 \
-n2 -m20 -s40 -g2000 -2K50m --secondary=no /path/to/reference/GRChg38.fq fq/${ID}_rename.fq  > bam/${ID}_minimap.sam

python retag.py --sam bam/${ID}_minimap.sam --out bam/${ID}_minitagged.bam
samtools sort -@ 10 bam/${ID}_minitagged.bam bam/${ID}_minitagged_sorted.bam
samtools index bam/${ID}_minitagged_sorted.bam
rm bam/${ID}_minimap.sam bam/${ID}_minitagged.bam
freebayes -f /path/to/reference/GRChg38.fq -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --limit-coverage 100000 bam/${ID}_minitagged.bam > vcf/${ID}.vcf
vartrix --umi --mapq 30 -b bam/${ID}_minitagged.bam -c ${DIR}/${ID}/outs/barcodes.tsv --scoring-method coverage --threads 16 --ref-matrix cts/${ID}_ref.mtx \
--out-matrix cts/${ID}_alt.mtx -v vcf/${ID}.vcf --fasta fq/${ID}_rename.fq
souporcell -a cts/${ID}_alt.mtx -r cts/${ID}_ref.mtx -b ${DIR}/${ID}/outs/barcodes.tsv -k 2 -t 8 > ${ID}_clusters.tsv
