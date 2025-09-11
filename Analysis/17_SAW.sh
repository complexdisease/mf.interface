ID=$1
image=$2
outdir=$3
cd ${outdir}
~/tools/saw-8.1.2/bin/saw count \
    --id=${ID} \
    --sn=${ID} \
    --omics=transcriptomics \
    --kit-version="Stereo-seq T FF V1.3" \
    --sequencing-type="PE75_50+100" \
    --chip-mask=/path/to/chip_mask/${ID}.barcodeToPos.h5 \
    --organism=human \
    --fastqs=/path/to/saw_fastq/${ID}/ \
    --reference=~/tools/saw-8.1.2/references/v8/human_transcriptome/\
    --image-tar=/path/to/saw_image/${image}
