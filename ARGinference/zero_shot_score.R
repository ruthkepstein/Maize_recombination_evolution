#extracting maize-specific hotspots
allchr_maizehotspots = allchr_whotspots %>% filter(hotspots == 'maize-specific')
#extracting teosinte-specific hotspots
allchr_teohotspots = allchr_whotspots %>% filter(hotspots == 'teosinte-specific')
#extracting shared hotspots
allchr_shared = allchr_whotspots %>% filter(hotspots == 'shared')
#extracting everything else
all_chr_all = allchr_whotspots %>% filter(hotspots == 'none')

#zero-shot score analysis, filtered by SNPs used in ARG-inference
zero_shot_score = read.table("/home2/rke27/Maize_Caduceus_zero_shot_hmp3filter_5cols_v4_sorted.bed")
colnames(zero_shot_score) = c("Chr", "Start", "End", "p", "Score")
#taking the bottom 5% of scores (most deleterious)
top_del = zero_shot_score %>% top_frac(-.05)

#reading in CDS sequences from B73v4 GFF
cds_bed = read.table("B73_v4_CDS_fixed.bed")
colnames(cds_bed)= c("Chr", "Start", "End")

#function to calculate zero-shot score per 1Mb interval
hotspot_score <- function(allchr, top_del, cds_bed) {
  
  results <- list()
  
  for(chr in unique(allchr$Chr)) {
    win_chr <- allchr[allchr$Chr == chr, ]
    top_chr <- top_del[top_del$Chr == chr, ]
    cds_chr <- cds_bed[cds_bed$Chr == chr, ]
    
    top_counts <- sapply(1:nrow(win_chr), function(i) {
      sum(top_chr$Start < win_chr$End[i] & top_chr$End > win_chr$Start[i])
    })
    
    #calculating median zero-shot score per interval
    median_scores <- sapply(1:nrow(win_chr), function(i) {
      overlapping <- top_chr[top_chr$Start < win_chr$End[i] & top_chr$End > win_chr$Start[i], ]
      if(nrow(overlapping) == 0) return(NA_real_)
      mean(overlapping$Score, na.rm = TRUE)
    })
    
    #calculating overlap of CDS sequences with each interval
    cds_bp <- sapply(1:nrow(win_chr), function(i) {
      overlaps <- cds_chr[cds_chr$Start < win_chr$End[i] & cds_chr$End > win_chr$Start[i], ]
      if(nrow(overlaps) == 0) return(0)
      overlap_lengths <- pmin(as.numeric(overlaps$End), as.numeric(win_chr$End[i])) -
        pmax(as.numeric(overlaps$Start), as.numeric(win_chr$Start[i]))
      overlap_lengths <- overlap_lengths[overlap_lengths > 0]
      if(length(overlap_lengths) == 0) return(0)
      sum(overlap_lengths, na.rm = TRUE)
    })
    
    #creating chromosome-specific df
    chr_res <- data.frame(
      Chr = win_chr$Chr,
      Start = win_chr$Start,
      End = win_chr$End,
      top_del = top_counts,
      median_score = median_scores,
      total_cds_bp = cds_bp
    )
    
    results[[as.character(chr)]] <- chr_res
  }
  
  #combining chromosomes
  final <- do.call(rbind, results)
  
  #calculating # of codons per interval by dividing # of bp by 3
  final$codons <- floor(final$total_cds_bp / 3)
  #load per codon
  final$load_per_codon <- ifelse(final$codons > 0, final$top_del / final$codons, NA)
  #score per codon
  final$score_per_codon <- ifelse(final$codons > 0, final$median_score / final$codons, NA)
  
  return(final)
}

#outputting scores in df
maize_score <- hotspot_score(allchr_maizehotspots, top_del, cds_bed)
teo_score <- hotspot_score(allchr_teohotspots, top_del, cds_bed)
shared_score <- hotspot_score(allchr_shared, top_del, cds_bed)
all_score <- hotspot_score(all_chr_all, top_del, cds_bed)

#t-test on each comparison
t.test(maize_score$score_per_codon, teo_score$score_per_codon)
t.test(shared_score$score_per_codon, teo_score$score_per_codon)
t.test(shared_score$score_per_codon, maize_score$score_per_codon)
