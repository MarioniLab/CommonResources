# Generating subread indices

## Common index combinations

Mouse mm10 genome by itself:

```sh
subread-buildindex -o mm10 /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa
```

Mouse mm10 genome with ERCC spike-ins:

```sh
subread-buildindex -o mm10_ERCC /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa ../sequences/spikes/ERCC92.fa
```

Human hg38 genome by itself:

```sh
subread-buildindex -o hg38 /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa
```

Human hg38 genome with ERCC spike-ins:

```sh
subread-buildindex -o hg38_ERCC /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa ../sequences/spikes/ERCC92.fa
```

## Genome builds for the spike-in project with Fernando Calero

This combines mm10 with both the ERCC and SIRV spike-ins:

```sh
subread-buildindex -o mm10_ERCC_SIRV /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
        ../sequences/spikes/ERCC92.fa ../sequences/spikes/SIRV_150601a.fasta 
```

... and an extra induced oncogene:

```sh
subread-buildindex -o mm10_ERCC_SIRV_onco /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
	../sequences/spikes/ERCC92.fa ../sequences/spikes/SIRV_150601a.fasta ../sequences/misc/CBFB-MYH11-mcherry.fa
```

Also building indices for the spike-ins by themselves, to assess specificity:

```sh
subread-buildindex -o ERCC ../sequences/spikes/ERCC92.fa 
subread-buildindex -o SIRV ../sequences/spikes/SIRV_150601a.fasta
```

See `/lustre/jmlab/lun01/SpikeIns` for more details.

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

## Genome builds for Lovorka Stojic's lncRNA project

This combines the hg38 genome with a plasmid sequence expressing dCas9-KRAB.

```sh
subread-buildindex -o hg38_cas9 ../sequences/misc/pHR-SFFV-dCas9-BFP-KRAB.fa \
    /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa
```

See `/lustre/jmlab/lun01/Odom/lncRNA_mitosis` for more details.

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

