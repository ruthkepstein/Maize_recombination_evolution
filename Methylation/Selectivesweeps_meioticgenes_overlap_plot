library(ggplot)

##plotting selective sweeps & meiotic genes in sweeps
library(GenomicRanges)
#lengths of maize chromosomes
chr_lengths <- data.frame(
  chromosome = paste0("chr", 1:10),
  length = c(307074717, 244480000, 235717834, 247040300, 223947240, 174083170, 182431542, 181172637, 160000002, 151082314)
)

genes <- data.frame(
  chromosome = c("chr6", "chr3", "chr4", "chr4", "chr4", "chr2", "chr3", "chr5", "chr4"),
  start = c(170757706, 194439248, 93519946, 104299642, 132079933, 204150197, 96131941, 60851572, 34374982),
  end = c(170774888, 194451193, 93539498, 104384163, 132107695, 204177280, 96172832, 60860717, 34377834),
  label = c("Rec8", "Exo1A", "Prd2", "Rad50", "Pms1", "Mlh3", "Figl1", "Mtop6b1", "Spo11-2")
)

sweeps = read.table("/home2/rke27/domestication_improvement_loci/domestication_Xu2020_v4_Teosinte_Maize_sorted_names.bed")
colnames(sweeps) = c("chromosome", "start", "end", "label")
sweeps = sweeps[-1,]

sweeps$type <- 'sweep'
genes$type <- 'gene'

gr_chr_lengths <- GRanges(
  seqnames = Rle(chr_lengths$chromosome),
  ranges = IRanges(start = 1, end = chr_lengths$length)
)

genes_gr <- GRanges(
  seqnames = genes$chromosome,
  ranges = IRanges(start = genes$start, end = genes$end),
  label = genes$label,
  type = genes$type
)

sweeps_gr <- GRanges(
  seqnames = sweeps$chromosome,
  ranges = IRanges(start = sweeps$start, end = sweeps$end),
  label = sweeps$label,
  type = sweeps$type
)

overlaps <- findOverlaps(sweeps_gr, genes_gr)
overlapping_regions <- pintersect(sweeps_gr[queryHits(overlaps)], genes_gr[subjectHits(overlaps)])

overlapping_df <- as.data.frame(overlapping_regions)

#ordering by chromosome # first
sweeps$chromosome <- factor(sweeps$chromosome,
                                      levels = paste0('chr', sort(as.numeric(gsub('chr', '', unique(sweeps$chromosome))), decreasing = TRUE)))

genes$chromosome <- factor(genes$chromosome,
                           levels = paste0('chr', sort(as.numeric(gsub('chr', '', unique(genes$chromosome))), decreasing = TRUE)))

overlapping_df$seqnames <- factor(overlapping_df$seqnames,
                                  levels = paste0('chr', sort(as.numeric(gsub('chr', '', unique(overlapping_df$seqnames))), decreasing = TRUE))
)

p <- ggplot() +
  #selective sweeps
  geom_segment(data = sweeps, aes(x = start/1000000, xend = end/1000000, y = chromosome, yend = chromosome),
               color = 'blue', size = 2, alpha = 0.4) +
  #genes
  geom_segment(data = genes, aes(x = start/1000000, xend = end/1000000, y = chromosome, yend = chromosome),
               color = 'green', size = 2, alpha = 0.4) +
  #overlapping regions
  geom_segment(data = overlapping_df, aes(x = (start)/1000000, xend = (end)/1000000, y = seqnames, yend = seqnames),
               color = 'red', size = 4, alpha = 1) +  # Thicker and more opaque
  #annotations with gene names for overlaps
  geom_text_repel(data = genes, aes(x = (start/1000000 + end/1000000) / 2, y = chromosome, label = label),
                  vjust = -0.8, color = 'red', fontface = 'bold', size = 3, nudge_y = 0.2, max.overlaps = Inf) +
  labs(title = "Domestication Selective Sweeps & Recombination-modifiers", x = "Position (Mb)", y = "Chromosome") +
  theme_minimal()

print(p)
