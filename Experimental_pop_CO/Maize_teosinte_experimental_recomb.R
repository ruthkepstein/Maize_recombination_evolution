Hi there, this is the final version of the teosinte project code all in one doc
Any questions can be sent to: jjwheeler4120@gmail.com or rke27@cornell.edu

This document has:
the landscape plotting with crossover positions

# Load required packages
library(ggplot2)
library(vroom)
library(scales)
library(dplyr)
library(ggvenn)
library(readxl)

m_intervals <- read.table("NAM_US_COs_v4_wind.bed", head = T)
colnames(m_intervals) = c("Chr", "Start", "End", "Extend","FullSampleName","length")
#removing NAs
m_intervals <- m_intervals[!is.na(m_intervals$Chr),]
m_intervals$length = m_intervals$End - m_intervals$Start
# TEOSINTE
t_intervals <- vroom("xo_imputed_Panzea_Teosinte_C2_and_PanzeaC2_Teo_Parents_withWI_sometaxa_removed.txt")

# MAIZE, only midpoint needed
m_intervals$mid = (m_intervals$End+m_intervals$Start)/2
# TEOSINTE
t_intervals$length = t_intervals$end-t_intervals$start
t_intervals$mid = (t_intervals$end+t_intervals$start)/2

# MAIZE
# this adds a count of the number of rows/crossovers with the same individual ID
m_ndup <- m_intervals %>% add_count(FullSampleName, name = "ndup")
m_filtered <- m_ndup[which(m_ndup$length <= 100000),]
m_ndup_filt = m_filtered[which(m_filtered$ndup <= 40),]

# TEOSINTE
# this adds a count of the number of rows/crossovers with the same individual ID
t_ndup <- t_intervals %>% add_count(taxon, name = "ndup")
t_ndup_filtered = t_ndup[which(t_ndup$length <= 100000),]
t_ndup_filt = t_ndup_filtered[which(t_ndup_filtered$ndup <= 30),]

# MAIZE
m_chr1 <- m_filtered[which(m_filtered$Chr==1),]
m_chr2 <- m_filtered[which(m_filtered$Chr==2),]
m_chr3 <- m_filtered[which(m_filtered$Chr==3),]
m_chr4 <- m_filtered[which(m_filtered$Chr==4),]
m_chr5 <- m_filtered[which(m_filtered$Chr==5),]
m_chr6 <- m_filtered[which(m_filtered$Chr==6),]
m_chr7 <- m_filtered[which(m_filtered$Chr==7),]
m_chr8 <- m_filtered[which(m_filtered$Chr==8),]
m_chr9 <- m_filtered[which(m_filtered$Chr==9),]
m_chr10 <- m_filtered[which(m_filtered$Chr==10),]

# TEOSINTE
t_chr1 <- t_ndup_filtered[which(t_ndup_filtered$chr==1),]
t_chr2 <- t_ndup_filtered[which(t_ndup_filtered$chr==2),]
t_chr3 <- t_ndup_filtered[which(t_ndup_filtered$chr==3),]
t_chr4 <- t_ndup_filtered[which(t_ndup_filtered$chr==4),]
t_chr5 <- t_ndup_filtered[which(t_ndup_filtered$chr==5),]
t_chr6 <- t_ndup_filtered[which(t_ndup_filtered$chr==6),]
t_chr7 <- t_ndup_filtered[which(t_ndup_filtered$chr==7),]
t_chr8 <- t_ndup_filtered[which(t_ndup_filtered$chr==8),]
t_chr9 <- t_ndup_filtered[which(t_ndup_filtered$chr==9),]
t_chr10 <- t_ndup_filtered[which(t_ndup_filtered$chr==10),]


# Plot Maize data - crossover count
# you can change bin size by altering "binwidth" in geom_histogram()
# CHROMOSOME 1
m_chr1_plot <- ggplot(m_chr1, aes(x = mid)) + # x-axis = midpoint position
  geom_histogram(binwidth = 1000000, # interval size = 60 kb
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 1 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw() # remove ugly ass gray background
m_chr1_plot
# Save the plot data
m_chr1_data <- ggplot_build(m_chr1_plot)$data[[1]]
# Calculate recombination rate in each interval
pop_size <- nrow(m_chr1[!duplicated(m_chr1$FullSampleName),]) # record population size
for(i in 1:nrow(m_chr1_data)){
  m_chr1_data[i,19] <- (((m_chr1_data[i,2] / pop_size)*100)/1)
} # Rate = ((#COs in interval / pop. size) * 100) / interval length in Mb)
colnames(m_chr1_data)[19] <- "r_rate"

# CHROMOSOME 2
m_chr2_plot <- ggplot(m_chr2, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 2 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr2_plot
m_chr2_data <- ggplot_build(m_chr2_plot)$data[[1]]
pop_size <- nrow(m_chr2[!duplicated(m_chr2$FullSampleName),])
for(i in 1:nrow(m_chr2_data)){
  m_chr2_data[i,19] <- (((m_chr2_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr2_data)[19] <- "r_rate"

# CHROMOSOME 3
m_chr3_plot <- ggplot(m_chr3, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 3 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr3_plot
m_chr3_data <- ggplot_build(m_chr3_plot)$data[[1]]
pop_size <- nrow(m_chr3[!duplicated(m_chr3$FullSampleName),])
for(i in 1:nrow(m_chr3_data)){
  m_chr3_data[i,19] <- (((m_chr3_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr3_data)[19] <- "r_rate"

# CHROMOSOME 4
m_chr4_plot <- ggplot(m_chr4, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 4 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr4_plot
m_chr4_data <- ggplot_build(m_chr4_plot)$data[[1]]
pop_size <- nrow(m_chr4[!duplicated(m_chr4$FullSampleName),])
for(i in 1:nrow(m_chr4_data)){
  m_chr4_data[i,19] <- (((m_chr4_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr4_data)[19] <- "r_rate"

# CHROMOSOME 5
m_chr5_plot <- ggplot(m_chr5, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 5 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr5_plot
m_chr5_data <- ggplot_build(m_chr5_plot)$data[[1]]
pop_size <- nrow(m_chr5[!duplicated(m_chr5$FullSampleName),])
for(i in 1:nrow(m_chr5_data)){
  m_chr5_data[i,19] <- (((m_chr5_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr5_data)[19] <- "r_rate"

# CHROMOSOME 6
m_chr6_plot <- ggplot(m_chr6, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 6 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr6_plot
m_chr6_data <- ggplot_build(m_chr6_plot)$data[[1]]
pop_size <- nrow(m_chr6[!duplicated(m_chr6$FullSampleName),])
for(i in 1:nrow(m_chr6_data)){
  m_chr6_data[i,19] <- (((m_chr6_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr6_data)[19] <- "r_rate"

# CHROMOSOME 7
m_chr7_plot <- ggplot(m_chr7, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 7 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr7_plot
m_chr7_data <- ggplot_build(m_chr7_plot)$data[[1]]
pop_size <- nrow(m_chr7[!duplicated(m_chr7$FullSampleName),])
for(i in 1:nrow(m_chr7_data)){
  m_chr7_data[i,19] <- (((m_chr7_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr7_data)[19] <- "r_rate"

# CHROMOSOME 8
m_chr8_plot <- ggplot(m_chr8, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 8 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr8_plot
m_chr8_data <- ggplot_build(m_chr8_plot)$data[[1]]
pop_size <- nrow(m_chr8[!duplicated(m_chr8$FullSampleName),])
for(i in 1:nrow(m_chr8_data)){
  m_chr8_data[i,19] <- (((m_chr8_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr8_data)[19] <- "r_rate"

# CHROMOSOME 9
m_chr9_plot <- ggplot(m_chr9, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 9 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr9_plot
m_chr9_data <- ggplot_build(m_chr9_plot)$data[[1]]
pop_size <- nrow(m_chr9[!duplicated(m_chr9$FullSampleName),])
for(i in 1:nrow(m_chr9_data)){
  m_chr9_data[i,19] <- (((m_chr9_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr9_data)[19] <- "r_rate"

# CHROMOSOME 10
m_chr10_plot <- ggplot(m_chr10, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkred", fill = "darkred", alpha = 0.5) +
  ggtitle("Maize Chromosome 10 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
m_chr10_plot
m_chr10_data <- ggplot_build(m_chr10_plot)$data[[1]]
pop_size <- nrow(m_chr10[!duplicated(m_chr10$FullSampleName),])
for(i in 1:nrow(m_chr10_data)){
  m_chr10_data[i,19] <- (((m_chr10_data[i,2] / pop_size) * 100)/1)
}
colnames(m_chr10_data)[19] <- "r_rate"

# Plot Teosinte data - crossover count

# CHROMOSOME 1
t_chr1_plot <- ggplot(t_chr1, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 1 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr1_plot
t_chr1_data <- ggplot_build(t_chr1_plot)$data[[1]]
pop_size <- nrow(t_chr1[!duplicated(t_chr1$taxon),])
for(i in 1:nrow(t_chr1_data)){
  t_chr1_data[i,19] <- (((t_chr1_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr1_data)[19] <- "r_rate"

# CHROMOSOME 2
t_chr2_plot <- ggplot(t_chr2, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 2 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr2_plot
t_chr2_data <- ggplot_build(t_chr2_plot)$data[[1]]
pop_size <- nrow(t_chr2[!duplicated(t_chr2$taxon),])
for(i in 1:nrow(t_chr2_data)){
  t_chr2_data[i,19] <- (((t_chr2_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr2_data)[19] <- "r_rate"

# CHROMOSOME 3
t_chr3_plot <- ggplot(t_chr3, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 3 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr3_plot
t_chr3_data <- ggplot_build(t_chr3_plot)$data[[1]]
pop_size <- nrow(t_chr3[!duplicated(t_chr3$taxon),])
for(i in 1:nrow(t_chr3_data)){
  t_chr3_data[i,19] <- (((t_chr3_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr3_data)[19] <- "r_rate"

# CHROMOSOME 4
t_chr4_plot <- ggplot(t_chr4, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 4 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr4_plot
t_chr4_data <- ggplot_build(t_chr4_plot)$data[[1]]
pop_size <- nrow(t_chr4[!duplicated(t_chr4$taxon),])
for(i in 1:nrow(t_chr4_data)){
  t_chr4_data[i,19] <- (((t_chr4_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr4_data)[19] <- "r_rate"

# CHROMOSOME 5
t_chr5_plot <- ggplot(t_chr5, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 5 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr5_plot
t_chr5_data <- ggplot_build(t_chr5_plot)$data[[1]]
pop_size <- nrow(t_chr5[!duplicated(t_chr5$taxon),])
for(i in 1:nrow(t_chr5_data)){
  t_chr5_data[i,19] <- (((t_chr5_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr5_data)[19] <- "r_rate"

# CHROMOSOME 6
t_chr6_plot <- ggplot(t_chr6, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 6 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr6_plot
t_chr6_data <- ggplot_build(t_chr6_plot)$data[[1]]
pop_size <- nrow(t_chr6[!duplicated(t_chr6$taxon),])
for(i in 1:nrow(t_chr6_data)){
  t_chr6_data[i,19] <- (((t_chr6_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr6_data)[19] <- "r_rate"

# CHROMOSOME 7
t_chr7_plot <- ggplot(t_chr7, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 7 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr7_plot
t_chr7_data <- ggplot_build(t_chr7_plot)$data[[1]]
pop_size <- nrow(t_chr7[!duplicated(t_chr7$taxon),])
for(i in 1:nrow(t_chr7_data)){
  t_chr7_data[i,19] <- (((t_chr7_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr7_data)[19] <- "r_rate"

# CHROMOSOME 8
t_chr8_plot <- ggplot(t_chr8, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 8 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr8_plot
t_chr8_data <- ggplot_build(t_chr8_plot)$data[[1]]
pop_size <- nrow(t_chr8[!duplicated(t_chr8$taxon),])
for(i in 1:nrow(t_chr8_data)){
  t_chr8_data[i,19] <- (((t_chr8_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr8_data)[19] <- "r_rate"

# CHROMOSOME 9
t_chr9_plot <- ggplot(t_chr9, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 9 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr9_plot
t_chr9_data <- ggplot_build(t_chr9_plot)$data[[1]]
pop_size <- nrow(t_chr9[!duplicated(t_chr9$taxon),])
for(i in 1:nrow(t_chr9_data)){
  t_chr9_data[i,19] <- (((t_chr9_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr9_data)[19] <- "r_rate"

# CHROMOSOME 10
t_chr10_plot <- ggplot(t_chr10, aes(x = mid)) +
  geom_histogram(binwidth = 1000000, 
                 color = "darkblue", fill = "darkblue", alpha = 0.5) +
  ggtitle("Teosinte Chromosome 10 Crossovers") +
  xlab("Position (bp)") + 
  theme_bw()
t_chr10_plot
t_chr10_data <- ggplot_build(t_chr10_plot)$data[[1]]
pop_size <- nrow(t_chr10[!duplicated(t_chr10$taxon),])
for(i in 1:nrow(t_chr10_data)){
  t_chr10_data[i,19] <- (((t_chr10_data[i,2] / pop_size) * 100)/1)
}
colnames(t_chr10_data)[19] <- "r_rate"
