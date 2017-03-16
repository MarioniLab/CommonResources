# Generating STAR indices

## Genome builds for Omer Ziv's Zika project

```{r}
mkdir temp
cat ../../annotation/processed/hg38.gtf ../../annotation/original/ZIKV_H.sapiens_Brazil_PE243_2015-1.gtf > temp/combined.gtf
mkdir hg38_zikv
bsub -R "rusage[mem=40000]" -n 4 -o log.out -e log.err STAR --runMode genomeGenerate --runThreadN 4 --genomeDir hg38_zikv \
    --genomeFastaFiles /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa ../sequences/misc/ZIKV_H.sapiens_Brazil_PE243_2015-1.fa \
    --sjdbGTFfile temp/combined.gtf --sjdbOverhang 100
```
