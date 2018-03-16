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

## Generating a _Xenopus laevis_ package for Elia Benito-Gutierrez's project

Building the _X. laevis_ genome, using a processed GTF file.

```sh
bsub -R "rusage[mem=32000]" -n 4 -e "xla.err" -o "xla.out" \
    /home/mib-cri/software/10Xgenomics/cellranger-1.2.0/cellranger mkref \
    --nthreads=4 --genome=xla \
    --fasta=/lustre/reference_data/mib-cri/reference_genomes/xenopus_laevis/JGI_9.1/fasta/xla.JGI_9.1.fa \
    --genes=/lustre/jmlab/resources/annotation/processed/XL_9.1_v1.8.3.2.gtf
```

Building the _X. tropicalis_ genome, to use instead of _laevis_ (trading off some bias for a more reliable build).

```sh
bsub -R "rusage[mem=32000]" -n 4 -e "xtr.err" -o "xtr.out" \
    /home/mib-cri/software/10Xgenomics/cellranger-1.2.0/cellranger mkref \
    --nthreads=4 --genome=xtr \
    --fasta=/lustre/reference_data/mib-cri/reference_genomes/xenopus_tropicalis/JGI_4.2/fasta/xtr.JGI_4.2.fa \
    --genes=/lustre/jmlab/resources/annotation/processed/xtr4.2.gtf
```

## Generating an Amphioxus package for Elia Benito-Gutierrez

```sh
#!/bin/bash
#SBATCH -n 4 
#SBATCH --mem 16000 
#SBATCH -e blan.err 
#SBATCH -o blan.out
/Users/bioinformatics/software/cellranger/cellranger-2.1.0/cellranger mkref \
    --nthreads=4 --genome=blan \
    --fasta=/scratchb/jmlab/resources/genomes/sequences/B.LAN_REFERENCE/Bl71nemr.fa \
    --genes=/scratchb/jmlab/resources/annotation/processed/B.LAN_exons.gtf
```


