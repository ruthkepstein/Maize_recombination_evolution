library(tidyverse)
library(stringr)
library(dplyr)
recombrate_maize_teo = function(r,m){
  ##where r is chr # and m is resolution in bp
  #this function reads in recomb & branchlen file from ARGweaver
  #returns recombination rates per m resolution for teosinte and maize
  setwd(paste("/home2/rke27/argweaver_chr",r,"_60maizelines_phased",sep =""))
  maize_arg = read.table(paste("maize51_282set_v4_chr",r,"_combined_recombrate_perMCMC_recent.bed", sep =""), fill = TRUE, row.names = NULL, header = F)
  colnames(maize_arg) = c("Chr", "Start", "End", "MCMCiter", "rate")
  maize_arg_sorted <- maize_arg[order(maize_arg$Start),]
  
  setwd(paste("/home2/rke27/argweaver_chr",r,"_teosintes_phased",sep =""))
  teo_arg = read.table(paste("teosintes_v4_chr",r,"_combined_recombrate_perMCMC_recent.bed", sep =""), fill = TRUE, row.names = NULL, header = F)
  colnames(teo_arg) = c("Chr", "Start", "End", "MCMCiter", "rate")
  teo_arg_sorted <- teo_arg[order(teo_arg$Start),]
  
  maize_arg_sorted = na.omit(maize_arg_sorted)
  teo_arg_sorted = na.omit(teo_arg_sorted)
  teo_arg_sorted <- subset(teo_arg_sorted, is.finite(rate))
  maize_arg_sorted <- subset(maize_arg_sorted, is.finite(rate))
  
  #now we must normalize by bp to get rate per 1bp
  maize_arg_sorted$length = maize_arg_sorted$End - maize_arg_sorted$Start
  teo_arg_sorted$length = teo_arg_sorted$End - teo_arg_sorted$Start
  maize_arg_sorted$rate_bp = maize_arg_sorted$rate/maize_arg_sorted$length
  teo_arg_sorted$rate_bp = teo_arg_sorted$rate/teo_arg_sorted$length
  
  breaks = seq(min(maize_arg_sorted$Start), max(maize_arg_sorted$Start), by = m)
  maize2 = maize_arg_sorted %>% mutate(bins = cut(Start, breaks = breaks))%>% 
    group_by(bins) %>% summarise(mean_rate = mean(rate_bp))
  teo2 = teo_arg_sorted %>% mutate(bins = cut(Start, breaks = breaks))%>% 
    group_by(bins) %>% summarise(mean_rate = mean(rate_bp))
  
  #combined = maize2
  combined = left_join(maize2, teo2, by = "bins")
  colnames(combined)=c("bin1", "maize_arg", "teo_arg")
  combined = tidyr::separate(combined,bin1,into = c("Start","End"),sep = ",",remove = FALSE,extra = "merge")
  combined$Start = str_sub(combined$Start,2,-1)
  combined$End = gsub('.{1}$','', combined$End)
  combined$Start = as.numeric(combined$Start)
  combined$End = as.numeric(combined$End)
  combined$Chr = r
  
  options(scipen = 3)
  write.table(combined[,c(6,2,3,4,5)], file = paste("/home2/rke27/crossover_datasets/Maize_Teosinte_ARGinferred_v4_chr",r,"_",m,"bp.bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  return(combined)
}

chr1_combined = recombrate_maize_teo(1, 1000000)
chr2_combined = recombrate_maize_teo(2, 1000000)
chr3_combined = recombrate_maize_teo(3, 1000000)
chr4_combined = recombrate_maize_teo(4, 1000000)
chr5_combined = recombrate_maize_teo(5, 1000000)
chr6_combined = recombrate_maize_teo(6, 1000000)
chr7_combined = recombrate_maize_teo(7, 1000000)
chr8_combined = recombrate_maize_teo(8, 1000000)
chr9_combined = recombrate_maize_teo(9, 1000000)
chr10_combined = recombrate_maize_teo(10, 1000000)

allchrs = rbind(chr1_combined, chr2_combined, chr3_combined, chr4_combined,
                chr5_combined, chr6_combined, chr7_combined,
                chr8_combined, chr9_combined, chr10_combined)

#comparing the maize and teosinte recombination landscapes
library(IDPmisc)
library(ggplot2)
combo_func = function(combined, r, s, num, m){
  #r is the chromosome # and s is the sd away from mean you want (1, 2 or 3 sd)
  #num is the top fraction of recombining regions (0.1 = 10%)
  #m is the bin length in kb (m = 100 = 100,000bp)
  #lets first plot distribution of maize & teosinte recombination rates across chromosome
  print(ggplot(combined, aes(x = Start/1000000)) +
    geom_line(aes(y = maize_arg, color = '#CC6677'), alpha = 0.5) +
    geom_line(aes(y = teo_arg, color = '#332288'), alpha = 0.5) +
    labs(x = paste("Chromosome",r," (Mb)",sep=""),
         y = "Recombination Rate") +
    scale_fill_identity(guide = 'legend') + theme_minimal() + theme(
      legend.position = c(.85, .85),
      legend.box.just = "right") +
    scale_colour_manual(name = "Population", values=c('#CC6677'='#CC6677','#332288'='#332288'), labels = c('teosinte','maize')))
  #calculating difference between maize & teosinte per bin
  combined$change = (combined$maize_arg-combined$teo_arg)
  #plotting to see difference between intervals
  plot(combined$Start, combined$change, type = 'l',xlab = "Maize - Teosinte")
  combined = na.omit(combined)
  #is normally distributed
  mean(combined$change)
  sd(combined$change)
  #want #s sd or greater from mean
  lower = mean(combined$change)-(sd(combined$change)*s)
  upper = mean(combined$change)+(sd(combined$change)*s)
  
  hist(combined$change, breaks = 80,xlab = "Maize - Teosinte")
  abline(v=lower, col = '#CC6677')
  abline(v=upper, col = '#332288')
  
  combined$Change="No Change"
  # Set new column values to appropriate colours
  combined$Change[combined$change>=upper]="Maize DRR"
  combined$Change[combined$change<=lower]="Teosinte DRR"
  print(ggplot(combined, aes(x = Start/1000000, y = change, col = Change)) + 
          geom_point(size = 0.8) +  # By default, show.legend is TRUE, so you can omit it
          xlab('Chromosome 1 (Mb)') + 
          ylab('Difference in recombination rates') +
          scale_color_manual(values = c("#332288", "grey", "#CC6677")) + 
          theme_bw() +
          theme(
            legend.position = c(0.5, 0.2),
            legend.box.just = "right"
          )+labs(color = NULL))
  
  setwd('/home2/rke27/crossover_datasets')
  maize_high_recombining_regions = combined[which(combined$Change == 'Hotter in Maize'), ]
  maize_high_recombining_regions$Chr = r
  maize_DRRs = maize_high_recombining_regions %>% select(Chr, Start, End, maize_arg)
  write.table(maize_DRRs, file = paste("Maize_DRRs_chr",r,"_",m,"kb_",s,"sd.bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  teo_high_recombining_regions = combined[which(combined$Change == 'Hotter in Teosinte'), ]
  teo_high_recombining_regions$Chr = r
  teo_DRRs = teo_high_recombining_regions %>% select(Chr, Start, End, teo_arg)
  write.table(teo_DRRs, file = paste("Teo_DRRs_chr",r,"_",m,"kb_",s,"sd.bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  ##looking at teosinte hotspots, regardless of being a DRR
  teo_arglen_hot = combined %>% top_frac(num, teo_arg)
  teo_arglen_hot$Chr = r
  teo_hotspots = teo_arglen_hot %>% select(Chr, Start, End, teo_arg)
  write.table(teo_hotspots, file = paste("teosinte_hotspots_top",num*100,"%_chr",r, ".bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  ##looking at maize hotspots, regardless of being a DRR
  frac = round(num*nrow(combined))
  maize_arglen_cold = combined %>% top_n(-frac, maize_arg)
  maize_arglen_cold$Chr = r
  maize_coldspots = maize_arglen_cold %>% select(Chr, Start, End, maize_arg)
  write.table(maize_coldspots, file = paste("maize_coldspots_top",num*100,"%_chr",r, ".bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  ##looking at maize hotspots, regardless of being a DRR
  maize_arglen_hot = combined %>% top_frac(num, maize_arg)
  maize_arglen_hot$Chr = r
  maize_hotspots = maize_arglen_hot %>% select(Chr, Start, End, maize_arg)
  write.table(maize_hotspots, file = paste("maize_hotspots_top",num*100,"%_chr",r, ".bed", sep = ""), col.names = F, row.names = F, sep = '\t')
  
  ##plotting where hotspots are on chromosomes & identifying if they are shared between maize & teosinte
  combined$h <- ifelse(combined$Start %in% maize_hotspots$Start, 'maize-specific', 'none')
  combined$h2 <- ifelse(combined$Start %in% teo_hotspots$Start, 'teosinte-specific', 'none')
  combined2 = combined %>% mutate(hotspots = case_when(h == 'maize-specific' & h2 == 'teosinte-specific' ~ 'shared',
                                                      h == 'maize-specific' & h2 == 'none' ~ 'maize-specific',
                                                      h == 'none' & h2 == 'teosinte-specific' ~ 'teosinte-specific',
                                             h == 'none' & h2 == 'none' ~ 'none'))
  shared = combined2 %>% count(hotspots == 'shared')
  teo = combined2 %>% count(hotspots == 'teosinte-specific')
  maize = combined2 %>% count(hotspots == 'maize-specific')
  print(ggplot(combined2, aes(x = Start/1000000, y = maize_arg, col = hotspots)) + 
          geom_point(size = 0.8) +
          xlab(paste('Chromosome ', r, ' (Mb)', sep = '')) + 
          ylab('Recombination Rate') +
          scale_color_manual(values = c("#332288", 'grey', '#A020F0', '#CC6677')) + 
          theme_bw() +
          theme(
            legend.position = c(0.5, 0.8),
            legend.box.just = "right")+labs(color = NULL))
  return_list = list("#shared" = shared[[2]][2], "#teosinte-specific" = teo[[2]][2], "#maize-specific" = maize[[2]][2])
  percent_shared = shared[[2]][2]/(shared[[2]][2] + maize[[2]][2])
  return(combined2)
}

chr1_whotspots = combo_func(chr1_combined, r = 1, s = 2, num = 0.1, m = 1000)
chr2_whotspots = combo_func(chr2_combined, r = 2, s = 2, num = 0.1, m = 100)
chr3_whotspots = combo_func(chr3_combined, r = 3, s = 2, num = 0.1, m = 1000)
chr4_whotspots = combo_func(chr4_combined, r = 4, s = 2, num = 0.1, m = 1000)
chr5_whotspots = combo_func(chr5_combined, r = 5, s = 2, num = 0.1, m = 1000)
chr6_whotspots = combo_func(chr6_combined, r = 6, s = 2, num = 0.1, m = 1000)
chr7_whotspots = combo_func(chr7_combined, r = 7, s = 2, num = 0.1, m = 1000)
chr8_whotspots = combo_func(chr8_combined, r = 8, s = 2, num = 0.1, m = 1000)
chr9_whotspots = combo_func(chr9_combined, r = 9, s = 2, num = 0.1, m = 1000)
chr10_whotspots = combo_func(chr10_combined, r = 10, s = 2, num = 0.1, m = 1000)

allchr_whotspots = rbind(chr1_whotspots, chr3_whotspots, chr4_whotspots, chr5_whotspots,
                         chr6_whotspots, chr7_whotspots, chr8_whotspots, chr9_whotspots,
                         chr10_whotspots)

write.table(allchr_whotspots[,c(6,2,3,4)], file = "/home2/rke27/crossover_datasets/Maize_DRR_hotspots.bed", col.names = F, sep = "\t", row.names = F)
write.table(allchr_whotspots[,c(6,2,3,5)], file = "/home2/rke27/crossover_datasets/Teosinte_DRR_hotspots.bed", col.names = F, sep = "\t", row.names = F)
