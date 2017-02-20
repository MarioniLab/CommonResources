# This script takes the hg38 RepeatMasker calls and converts it to a GTF file.
# The calls themselves were derived from http://www.repeatmasker.org/species/hg.html.

system('cat hg38.fa.out | sed "s/^  *//g" | sed "s/  */,/g" > blah.out')
repeat.data <- read.csv("blah.out", skip=2, fill=TRUE, stringsAsFactors=FALSE, header=FALSE, 
                        colClasses=list(NULL, NULL, NULL, NULL, "character", "integer", "integer", NULL, NULL, "character", NULL, NULL, NULL, NULL, NULL))
unlink("blah.out")

require(GenomicRanges)
rep.ranges <- GRanges(repeat.data[,1], IRanges(repeat.data[,2], repeat.data[,3]))

#reduced <- reduce(rep.ranges, min.gapwidth=100)
#olap <- findOverlaps(rep.ranges, reduced, select="first")
#re.name <- lapply(split(repeat.data[,4], olap), FUN=paste0, collapse="|")
#names(reduced)[as.integer(names(re.name))] <- unlist(re.name)

require(rtracklayer)
rep.ranges$type <- "exon"
rep.ranges$gene_id <- paste0(repeat.data[,4], "|", seq_along(rep.ranges))
export(rep.ranges, format="gtf", con="repeats_hg38.gtf")


