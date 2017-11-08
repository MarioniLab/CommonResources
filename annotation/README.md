# Original sources of all annotation

__Type__ | __Source__
--- | ---
ERCC92  |   http://www.thermofisher.com/order/catalog/product/4456739
CBFB-MYH11-mcherry  |   manually constructed, based on information from Fernando Calero
Homo_sapiens.GRCh37.75.gtf.gz	|	http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
Homo_sapiens.GRCh38.83.gtf.gz	|	http://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/
Homo_sapiens.GRCh38.90.gtf.gz  	|	http://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/
Mus_musculus.GRCm38.82.gtf.gz	|	http://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/
Mus_musculus.GRCm38.90.gtf.gz	|	http://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/
SIRV_C_150601a  |	https://www.lexogen.com/sirvs/
cas9_pHR_approx |   manually constructed by examining Cas9-coding domain in https://www.addgene.org/46911/
repeats/hg38.fa.out |	http://www.repeatmasker.org/species/hg.html
B.LAN_REFERENCE/*.gtf	|	supplied by Elia Benito-Gutierrez, via the EBI servers (Dec 16, 2016)
XL_9.1_v1.8.3.2.allTranscripts.gff3.gz  |   ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/1.8.3.2/
Xenopus_tropicalis.JGI_4.2.90.gtf.gz    |   ftp://ftp.ensembl.org/pub/release-90/gtf/xenopus_tropicalis/Xenopus_tropicalis.JGI_4.2.90.gtf.gz
ZIKV_H.sapiens_Brazil_PE243_2015-1.gtf  |   manually constructed as containing the entire Zika genome

# Processing Ensembl annotation

Ensembl GTF's don't put `chr` in front of their chromosome names.
We put them back in to avoid problems when matching up with FASTA sequence names.
Similarly, the mitochondrial chromosome is named as `MT`, which needs to be changed to `M`.
Finally, we only keep regions annotated as exons to avoid pulling down the coding region, start/stop codons, etc.

```sh
zcat original/Mus_musculus.GRCm38.90.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > processed/Mus_musculus.GRCm38.90.gtf
```

We do the same for hg38.

```sh
zcat original/Homo_sapiens.GRCh38.90.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > processed/Homo_sapiens.GRCh38.90.gtf
```

You can check proper naming of the chromosomes with `cut -f1 | uniq -c`.

There is no need for chromosomal name conversion with _Xenopus tropicalis_, but we do need to keep exons and to unzip the file.

```{r}
zcat original/Xenopus_tropicalis.JGI_4.2.90.gtf.gz | \
    awk '$3 == "exon"' > processed/Xenopus_tropicalis.JGI_4.2.90.gtf
```

# Converting RepeatMasker output to GTF

```sh
Rscript original/repeats/rm2gtf.R
mv original/repeats/repeats_hg38.gtf processed/
```

# Obtaining ID to name mappings for lancelet genes

The ID:name mappings are only available from the GTF file.
We extract them for convenience, to avoid having to repeatedly load the entire file to memory during analysis.

```sh
host=original/B.LAN_REFERENCE/Bla_annot-FINAL_v4_names-phCuff.gtf 
cat ${host} | cut -f9 | sed -r "s/.*gene_id ([^;]+);.*/\1/" > ids.txt
cat ${host} | cut -f9 | sed -r "s/.*gene_name ([^;]+);.*/\1/" > names.txt
paste ids.txt names.txt | uniq > processed/B.LAN_annotation.txt
rm ids.txt names.txt
```

# Processing _Xenopus laevis_ annotations

We need to extract only the exons (and UTRs) for counting.
This is done in R because we need map names to gene IDs as well.

```{r}
library(rtracklayer)
original <- import("original/XL_9.1_v1.8.3.2.allTranscripts.gff3.gz")
stopifnot(all(lengths(original$ID)==1L))
genes <- original[original$type=="mRNA"]
stopifnot(all(lengths(genes$Parent)==1L))
stopifnot(all(lengths(genes$Name)==1L))

# Making the GTF from only the CDS and UTRs.
keep <- original$type=="CDS" | grepl("UTR", original$type)
stripped <- original[keep,]
stripped$type <- "exon"
stripped$exon_id <- stripped$ID
stripped$transcript_id <- stripped$Parent
stripped$Name <- stripped$longest <- stripped$ID <- stripped$Parent <- NULL

# Mapping transcript IDs onto gene IDs.
m <- match(as.character(stripped$transcript_id), as.character(genes$ID))
stopifnot(all(!is.na(m)))
stripped$gene_id <- as.character(genes$Parent)[m]
stripped$gene_name <- as.character(genes$Name)[m]
export(stripped, con="processed/XL_9.1_v1.8.3.2.gtf")

# Mapping gene IDs to gene names.
write.table(data.frame(unlist(genes$Parent), unlist(genes$Name)),
    file="processed/XL_annotation.txt", sep="\t", row.names=FALSE, col.names=FALSE)
```
