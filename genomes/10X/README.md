# Generating 10X annotation packages

## Common use cases

Human (hg38) and mouse (mm10) packages can be downloaded from the 10X website.
They should be hosted centrally at `/lustre/reference_data/mib-cri/10XGenomics`.

## Generating tc1 packages for Celia Martinez's project

Adding the human chr21 to the mm10 genome (sequence and GTF) and calling `cellranger`.

```sh
cat /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
    /lustre/jmlab/resources/genomes/sequences/misc/hg38.chr21.fa > tc1.fa
cat /lustre/jmlab/resources/annotation/processed/hg38.gtf | grep "^chr21" | \
    cat /lustre/jmlab/resources/annotation/processed/mm10.gtf - > tc1.gtf
bsub -R "rusage[mem=32000]" -n 4 -e "tc1.err" -o "tc1.out" \
    /home/mib-cri/software/10Xgenomics/cellranger-1.2.0/cellranger mkref \
    --nthreads=4 --genome=tc1 --fasta=tc1.fa --genes=tc1.gtf 
```
