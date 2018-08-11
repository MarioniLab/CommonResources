# Generating bowtie2 indices

Human hg38 genome by itself:

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.hg38.out
#SBATCH -e log.hg38.err
#SBATCH -n 1
#SBATCH --mem 16000
/Users/bioinformatics/software/bowtie/bowtie2-2.3.3.1/bowtie2-build /mnt/scratcha/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa hg38
EOT
```

Lancelet draft genome with ERCC spike-ins.

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.blan.out
#SBATCH -e log.blan.err
#SBATCH -n 1
#SBATCH --mem 16000
/Users/bioinformatics/software/bowtie/bowtie2-2.3.3.1/bowtie2-build \
    /mnt/scratchb/jmlab/resources/genomes/sequences/B.LAN_REFERENCE/Bl71nemr.fa \
    /mnt/scratchb/jmlab/resources/genomes/sequences/spikes/ERCC92.fa -o BL71_ERCC 
EOT
```
