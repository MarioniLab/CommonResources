# Original sources of all annotation

__Type__ | __Source__
--- | ---
ERCC92  |   http://www.thermofisher.com/order/catalog/product/4456739
CBFB-MYH11-mcherry  |   manually constructed, based on information from Fernando Calero
Homo_sapiens.GRCh37.75.gtf  |	http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
Homo_sapiens.GRCh38.83.gtf  |	http://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/
Mus_musculus.GRCm38.82.gtf  |	http://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/
SIRV_C_150601a  |	https://www.lexogen.com/sirvs/
cas9_pHR_approx |   manually constructed by examining Cas9-coding domain in https://www.addgene.org/46911/
repeats/hg38.fa.out |	http://www.repeatmasker.org/species/hg.html

# Processing Ensembl annotation

Ensembl GTF's don't put `chr` in front of their chromosome names.
We put them back in to avoid problems when matching up with FASTA sequence names.
Similarly, the mitochondrial chromosome is named as `MT`, which needs to be changed to `M`.
Finally, we only keep regions annotated as exons to avoid pulling down the coding region, start/stop codons, etc.

```sh
cat original/Mus_musculus.GRCm38.82.gtf | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > processed/mm10.gtf
```

We do the same for hg38.

```sh
cat original/Homo_sapiens.GRCh38.83.gtf | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > processed/hg38.gtf
```

# Converting RepeatMasker output to GTF

```sh
Rscript original/repeats/rm2gtf.R
mv original/repeats/repeats_hg38.gtf processed/
```
