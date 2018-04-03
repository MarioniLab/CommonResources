# Generating subread indices

Mouse mm10 genome by itself:

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.mm10.out
#SBATCH -e log.mm10.err
#SBATCH -n 1
#SBATCH --mem 16000
subread-buildindex -o mm10 /scratchb/bioinformatics/reference_data/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa
EOT
```

Mouse mm10 genome with ERCC spike-ins:

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.mm10_ERCC.out
#SBATCH -e log.mm10_ERCC.err
#SBATCH -n 1
#SBATCH --mem 16000
subread-buildindex -o mm10_ERCC /scratchb/bioinformatics/reference_data/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa ../sequences/spikes/ERCC92.fa
EOT
```

Human hg38 genome by itself:

```sh
subread-buildindex -o hg38 /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa
```

Human hg38 genome with ERCC spike-ins:

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.hg38_ERCC.out
#SBATCH -e log.hg38_ERCC.err
#SBATCH -n 1
#SBATCH --mem 16000
subread-buildindex -o hg38_ERCC /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa ../sequences/spikes/ERCC92.fa
EOT
```

## Genome builds for Celia Martinez's tc1 project

This combines mm10 with chromosome 21 from hg38 and ERCC spike-ins

```sh
subread-buildindex -o mm10_human21_ERCC /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
	../sequences/misc/hg38.chr21.fa ../sequences/spikes/ERCC92.fa
```

See `/lustre/jmlab/lun01/Odom/trisomy21` for more details.

## Genome builds for Dan Greaves' HUSH project 

This constructs an index for repeat sequences.

```sh
subread-buildindex -o RepeatMasker ../sequences/repeats/RepeatMaskerLib.fasta
```

See `/lustre/jmlab/lun01/Lehner/HUSH_repeats` for more details.

## Genome builds for Elia Benito-Gutierrez's lancelet project.

This builds an index from the lancelet draft genome.

```sh
bsub -J builder -R "rusage[mem=10000]" \
    subread-buildindex ../B.LAN_REFERENCE/Bl71nemr.fa -o BL71
```

See `/lustre/jmlab/lun01/Benito-Gutierrez/devstage` for more details.

## Genome builds for Cut-and-run analysis

This builds an index from hg38 with spike-ins from yeast.

```sh
bsub -J builder -R "rusage[mem=10000]" \
    subread-buildindex -o hg38_yeast \
    /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa \
    /lustre/reference_data/mib-cri/reference_genomes/saccharomyces_cerevisiae/R64-1-1/fasta/sce.R64-1-1.fa
```

