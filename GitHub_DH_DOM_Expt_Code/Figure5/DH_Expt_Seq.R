library(tidyverse)
library(phyloseq)
library(readr)
library(dplyr)
library(ggplot2)


## Removed duplicate samples and Negs. Inverted DH56 (OSU moved samples so 1m was 24m etc) can also remove. Sample DH6 was labeled DH2-12 but that sample was not taken so relabled as DH2-5. Can also remove.
sample <- as.data.frame(read_csv('DH_Expt_Meta.csv'))
taxa <- as.data.frame(read_csv('DH_Expt_Taxa.csv'))
table <- as.data.frame(read_csv('DH_Expt_Table.csv'))
sample$Depth <- as.factor(sample$Depth)

## define the row names for the taxonomy table
row.names(table) <- table$Seq
row.names(taxa) <- taxa$Seq
row.names(sample) <- sample$Seq_ID

## remove the column Seq and ID
table <- table %>% select (-Seq)
taxa <- taxa %>% select (-Seq)
sample <- sample %>% select (-Seq_ID)

## make into matrices
table <- as.matrix(table)
taxa <- as.matrix(taxa)
sample <- as.data.frame(sample)

## transform to phyloseq objects
table <- otu_table(table, taxa_are_rows = TRUE)
taxa <- tax_table(taxa)
sample <- sample_data(sample)

## make the phyloseq bundle and save it
DH.expt <- phyloseq(table, taxa, sample)
saveRDS(DH.expt, file = "DH.expt.rds")

## ASV should be in more than four samples and its abundance should be greater than 0.0015 of the total
DH.expt = filter_taxa(DH.expt, function(x) sum(x > 3) > (0.0015*length(x)), TRUE)
DH.expt = prune_samples(sample_sums(DH.expt)>=1000, DH.expt)

## Plot the diversity indices for all the data, time-series = ts
dh.index1 <- plot_richness(DH.expt, x="Sample_ID", measures=c("Shannon", "Simpson"), color="Treatment")
dh.index1
## Remove the sample with low diversity E4
dh.expt <- prune_samples(sample_names(DH.expt) != 'DH36_E4', DH.expt)
##Replot diversity indices
dh.index2 <- plot_richness(dh.expt, x="Sample_ID", measures=c("Shannon", "Simpson"), color="Treatment")
dh.index2

## Plot NMDS for all the data - PC and Ominopore filters show differences - sequencing runs?
dh.nmds <- transform_sample_counts(dh.expt, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(dh.nmds, method="NMDS", distance="bray")
plot_ordination(dh.nmds, ord.nmds.bray, color="Treatment",  title="Bray NMDS")

dh.expt <- subset_samples(dh.expt, Type1 == 'Expt') ## subset to the samles involved in the experiment

## Plot bar charts for expt at Family level
profile_bar <- names(sort(taxa_sums(dh.expt), decreasing=TRUE))
profile_data <- transform_sample_counts(dh.expt, function(OTU) OTU/sum(OTU))
profile_plot <- prune_taxa(profile_bar, profile_data) 
bar_plot <- tax_glom(profile_plot, "Family")
plot_bar(bar_plot, x="Sample_ID", fill ="Family") + geom_bar(stat="identity")

## Export the bar plot to ggplot to allow formating
gg <- plot_bar(bar_plot, x="Sample_ID", fill ="Family")
ggdata1 <- gg$data

#Order the factors for Day
ggdata1 <- ggdata1 %>%
  arrange(Day2) %>%
  mutate(Day2 = factor(Day2, levels = c("Day 0","Day 2" , "Day 6","Day 9" ,"Day 12","Day 21")))

#Order the factors for Treatment
ggdata1 <- ggdata1 %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))
## reduce family to the more abundant ASVs using lin_focus_seq file 
##Order by Abundant Family
lin_focus_seq <- c("Actinomarinales","Anaerolineales","C. Kaiserbacteria","C. Moranbacteria","Chitinophagales","Chlorobiales","Dadabacteriales", "Desulfobacterales","Flavobacteriales","Ga0077536","Marine Group II", "Nitrosopumilales",
                   "Planctomycetales", "Pseudomonadales","Puniceispirillales","Rhodobacterales","Rhodospirillales","Rickettsiales","SAR11 clade","SAR202 clade", "SAR324 clade", "SAR406 clade", "Synechococcales","Other")

`%notin%` <- Negate(`%in%`)

ggdata1 <- ggdata1 %>%
  as_tibble() %>%
  arrange(Family) %>%
  relocate(Family) %>%
  mutate(Family = as.character(Family),
         Family = case_when(Family %notin% lin_focus_seq ~ "Other",
                            TRUE ~ Family),
         Family = factor(Family, levels = rev(lin_focus_seq)))

## Set colours
seq_colour <-c("grey70","green3","mediumpurple3","yellow" ,"maroon", "dodgerblue3", "tan","hotpink","lightpink", "darkseagreen1","coral", "deepskyblue", "forestgreen", "black", "mediumvioletred","gold","purple","white","red","cyan","wheat","brown4","gray40", "navy")

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command
seq_bar <- ggplot(data = ggdata1) + 
  geom_bar(aes(x = Day2, y = Abundance, fill = Family), stat = "Identity", position = "fill") +
  ylab("") + xlab("") +
  ## remove the grid background
  theme_bw() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Times New Roman")) +
  theme(axis.text.x = element_text(color = "black", family="Times New Roman", size = 18,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Times New Roman", size = 18)) +
  scale_fill_manual(values =seq_colour) + facet_grid(Filter ~Treatment) +
  theme(strip.text.x = element_text(color = "black", size = 20, face = "bold"))
seq_bar

dh.expt.nc <- subset_taxa(dh.expt, !Family == 'Chlorobiales') ## remove Chlorobiales
profile_bar <- names(sort(taxa_sums(dh.expt.nc), decreasing=TRUE)) [1:100]
profile_data <- transform_sample_counts(dh.expt.nc, function(OTU) OTU/sum(OTU))
profile_plot <- prune_taxa(profile_bar, profile_data) 
bar_plot <- tax_glom(profile_plot, "Family")
plot_bar(bar_plot, x="Sample_ID", fill ="Family") + geom_bar(stat="identity")

## Export the bar plot to ggplot to allow formating
gg.nc <- plot_bar(bar_plot, x="Sample_ID", fill ="Family")
ggdata.nc <- gg.nc$data

#Order the factors for Day
ggdata.nc <- ggdata.nc %>%
  arrange(Day2) %>%
  mutate(Day2 = factor(Day2, levels = c("Day 0","Day 2" , "Day 6","Day 9" ,"Day 12","Day 21")))

#Order the factors for Treatment
ggdata.nc <- ggdata.nc %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))
## reduce family to the more abundant ASVs using lin_focus_seq file 
##Order by Abundant Family
lin_focus_nc <- c("Actinomarinales","Anaerolineales","C. Kaiserbacteria","C. Moranbacteria","Chitinophagales","Dadabacteriales", "Desulfobacterales","Flavobacteriales","Ga0077536","Marine Group II", "Nitrosopumilales",
                   "Planctomycetales", "Pseudomonadales","Puniceispirillales","Rhodobacterales","Rhodospirillales","Rickettsiales","SAR11 clade","SAR202 clade", "SAR324 clade", "SAR406 clade", "Synechococcales","Other")

`%notin%` <- Negate(`%in%`)

ggdata.nc <- ggdata.nc %>%
  as_tibble() %>%
  arrange(Family) %>%
  relocate(Family) %>%
  mutate(Family = as.character(Family),
         Family = case_when(Family %notin% lin_focus_nc ~ "Other",
                            TRUE ~ Family),
         Family = factor(Family, levels = rev(lin_focus_nc)))

## Set colours
seq_colour_nc <-c("grey70","green3","mediumpurple3","yellow" ,"maroon", "dodgerblue3", "tan","hotpink","lightpink", "darkseagreen1","coral", "deepskyblue", "forestgreen", "black", "mediumvioletred","gold","purple","white","cyan","wheat","brown4","gray40", "navy")

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command
seq_bar_nc <- ggplot(data = ggdata.nc) + 
  geom_bar(aes(x = Day2, y = Abundance, fill = Family), stat = "Identity", position = "fill") +
  ylab("") + xlab("") +
  ## remove the grid background
  theme_bw() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Times New Roman")) +
  theme(axis.text.x = element_text(color = "black", family="Times New Roman", size = 18,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Times New Roman", size = 18)) +
  scale_fill_manual(values =seq_colour_nc) + facet_grid(Filter ~Treatment) +
  theme(strip.text.x = element_text(color = "black", size = 20, face = "bold"))
seq_bar_nc