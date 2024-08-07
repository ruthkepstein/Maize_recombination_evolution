#extracting maize-specific hotspots
allchr_maizehotspots = allchr_whotspots %>% filter(hotspots == 'maize-specific')
#extracting teosinte-specific hotspots
allchr_teohotspots = allchr_whotspots %>% filter(hotspots == 'teosinte-specific')
#zero-shot score analysis, filtered by SNPs used in ARG-inference
zero_shot_score = read.table("/home2/rke27/Maize_Caduceus_zero_shot_hmp3filter_5cols_v4_sorted.bed")
colnames(zero_shot_score) = c("Chr", "Start", "End", "p", "Score")

#function to calculate zero-shot score per 1Mb interval
hotspot_score = function(allchr){
  score <- allchr %>%
    inner_join(zero_shot_score, by = "Chr") %>%
    filter(Start.x < End.y & End.x > Start.y) %>%
    mutate(
      intersect_Start = pmax(Start.x, Start.y),
      intersect_End = pmin(End.x, End.y)
    ) %>%
    group_by(Chr, Start.x, End.x) %>%
    summarise(mean_Score = mean(Score), .groups = 'drop') %>%
    rename(Start = Start.x, End = End.x)
  return(score)
}

#this function is similar, but gives the distribution of zero-shot scores across hotspot intervals
all_means_score = function(allchr){
  intersected_means <- allchr %>%
    inner_join(zero_shot_score, by = "Chr") %>%
    filter(Start.x < End.y & End.x > Start.y) %>%
    mutate(
      Start = pmax(Start.x, Start.y),
      End = pmin(End.x, End.y)
    ) %>%
    group_by(Chr, Start, End) %>%
    summarise(mean_Score = mean(Score)) %>%
    ungroup() %>%
    select(Chr, Start, End, mean_Score)
  return(intersected_means)
}

maize_score = hotspot_score(allchr_maizehotspots)
teo_score = hotspot_score(allchr_teohotspots)

maize_score_dist = all_means_score(allchr_maizehotspots)
teo_score_dist = all_means_score(allchr_teohotspots)

#t-test to find significance between maize and teosinte high recombining regions
t.test(maize_score$mean_Score, teo_score$mean_Score)
t.test(maize_score_dist$mean_Score, teo_score_dist$mean_Score)
