DIR=$1
cd ${DIR}
for ID in */;
do
        mkdir -p ${DIR}/${ID}/copyscat
        python3 ~/tools/CopyscAT/process_fragment_file.py -i ${DIR}/${ID}/outs/atac_fragments.tsv.gz -o ${DIR}/${ID}/copyscat/processed_output.tsv -g /genome/copyscat_reference/hg38_chrom_sizes.tsv -b 1000000 -f 10000
done &
wait
