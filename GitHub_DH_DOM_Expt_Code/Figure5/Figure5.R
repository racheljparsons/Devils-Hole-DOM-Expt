library(tidyverse)
library(phyloseq)
library(readr)
library(dplyr)
library(ggplot2)
library("gridExtra")
library("ggpubr")

## Amino Acid Plot - Profile first
aa1 <- read.csv(file = "AA_profile.csv", header = TRUE)

#Order the factors for Day
aa1 <- aa1 %>%
  arrange(ID) %>%
  mutate(ID = factor(ID, levels = c("Surface", 
                                        "Deep")))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours
aa1_colour <-c("navy", "gold", "red4", "forestgreen", "purple4", "lightblue2", "lightgoldenrod1", "lightcoral", "lightgreen", "mediumpurple1", "orange", "mediumturquoise","lightpink","lightslategrey", "violetred", "black")

aa1_bar <- ggplot(data = aa1) +
  geom_bar(aes(x = ID, y = Conc, fill = Amino_Acid), stat = "Identity", position = "fill") +
  ylab("Molar Fraction") + xlab("") + ggtitle("A)") +
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
  theme(title = element_text(color = "black", family="Arial",size = 20)) +
  scale_fill_manual(values = aa1_colour) + facet_grid(~Profile) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 20, face = "bold"))
aa1_bar

# Amino Acid Plot - Expt
aa2 <- read.csv(file = "AA_plot2.csv", header = TRUE)

#Order the factors for Day
aa2 <- aa2 %>%
  arrange(Day2) %>%
  mutate(Day2 = factor(Day2, levels = c("Day 0", 
                                        "Day 2",
                                        "Day 6",
                                        "Day 12")))

#Order the factors for Treatment
aa2 <- aa2 %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours
aa2_colour <-c("navy", "gold", "red4", "forestgreen", "purple4", "lightblue2", "lightgoldenrod1", "lightcoral", "lightgreen", "mediumpurple1", "orange", "mediumturquoise","lightpink","lightslategrey", "violetred", "black")

aa2_bar <- ggplot(data = aa2) +
  geom_bar(aes(x = Day2, y = Conc, fill = Amino_Acid), stat = "Identity", position = "fill") +
  ylab("") + xlab("") +
  ## remove the grid background
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
  scale_fill_manual(values = aa2_colour) + facet_grid(~Treatment) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 20, face = "bold"))
aa2_bar

##Combining the two Amino Acid plots into one with one legend.
figure5a <- ggarrange(aa1_bar, aa2_bar, ncol = 2, nrow = 1,
                     common.legend = TRUE, legend = "right",
                     align = "hv", widths = c(1,3))
figure5a
ggsave("Figure5a.png", plot = figure5a, width = 12, height = 6, dpi = 300, units = "in")

## HR-DOM Plot - Profile
hrdom1 <- read.csv(file = "HR_profile.csv", header = TRUE)

#Order the factors for Depth
hrdom1 <- hrdom1 %>%
  arrange(Depth) %>%
  mutate(Depth = factor(Depth, levels = c("Surface", 
                                      "Deep")))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours Condensed Hydrocarbon	Lignin	Lipid	Protein1	Protein2	Carbohydrates	Protein Maya	Black Carbon	Polyphenols	Highly Unsaturated	Unsaturated Aliphatics	Peptides	Sugars	Saturated Fatty Acids	CRAM
hrdom1_colour <-c("black", "green3", "gold", "red" , "violet", "green4", "khaki", "blue", "orange", "purple3", "navy", "lightblue", "pink", "lightgreen", "tan")

hrdom1_bar <- ggplot(data = hrdom1) +
  geom_bar(aes(x = Depth, y = Peak_Area, fill = Compound), stat = "Identity", position = "fill") +
  ylab("Peak Area Fraction") + xlab("") +  ggtitle("B)") +
  ## remove the grid background
  ## remove the grid background
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
  theme(title = element_text(color = "black", family="Arial",size = 20)) +
  scale_fill_manual(values = hrdom1_colour) + facet_grid(~Type) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 20, face = "bold"))
hrdom1_bar

##HR-DOM Expt Plot
hrdom2 <- read.csv(file = "HR_expt.csv", header = TRUE)

#Order the factors for Day
hrdom2 <- hrdom2 %>%
  arrange(Day) %>%
  mutate(Day = factor(Day, levels = c("Day 0", 
                                      "Day 6",
                                      "Day 21")))

#Order the factors for Treatment
hrdom2 <- hrdom2 %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours Condensed Hydrocarbon	Lignin	Lipid	Protein1	Protein2	Carbohydrates	Protein Maya	Black Carbon	Polyphenols	Highly Unsaturated	Unsaturated Aliphatics	Peptides	Sugars	Saturated Fatty Acids	CRAM
hrdom2_colour <-c("black", "green3", "gold", "red" , "violet", "green4", "khaki", "blue", "orange", "purple3", "navy", "lightblue", "pink", "lightgreen", "tan")

hrdom2_bar <- ggplot(data = hrdom2) +
  geom_bar(aes(x = Day, y = Peak_Area, fill = Compound), stat = "Identity", position = "fill") +
  ylab("") + xlab("") +
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
  scale_fill_manual(values = hrdom2_colour) + facet_grid(~Treatment) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 20, face = "bold"))
hrdom2_bar

##Combining the two HR-DOM plots into one with one legend.
figure5b <- ggarrange(hrdom1_bar, hrdom2_bar, ncol = 2, nrow = 1,
                      common.legend = TRUE, legend = "right",
                      align = "hv", widths = c(1,3))
figure5b
ggsave("Figure5b.png", plot = figure5b, width = 12, height = 6, dpi = 300, units = "in")

## Create the sequencing bar plot
seq_profile <- readRDS("DH36.rds")
seq_profile = prune_samples(sample_sums(seq_profile)>=1000, seq_profile)


seq_profile <- subset_samples(seq_profile, Type == 'Expt') ## subset to profile
profile_bar <- names(sort(taxa_sums(seq_profile), decreasing=TRUE))
profile_data <- transform_sample_counts(seq_profile, function(OTU) OTU/sum(OTU))
profile_plot <- prune_taxa(profile_bar, profile_data) 
bar_plot <- tax_glom(profile_plot, "Family")
plot_bar(bar_plot, x="ID", fill ="Family") + geom_bar(stat="identity")

## Export the bar plot to ggplot to allow formatting
seq1 <- plot_bar(bar_plot, x="ID", fill ="Family")
seqdata1 <- seq1$data

#Order the factors for Day
seqdata1 <- seqdata1 %>%
  arrange(Depth2) %>%
  mutate(Depth2 = factor(Depth2, levels = c("Surface", 
                                            "Deep")))


##Order by Abundant Family
lin_focus_seq <- c("Actinomarinales","C. Kaiserbacteria","Chitinophagales","Chlorobiales","Dadabacteriales", "Desulfobacterales","Flavobacteriales","Ga0077536","Marine Group II", "Nitrosopumilales",
                    "Planctomycetales", "Pseudomonadales","Puniceispirillales","Rhodobacterales","Rhodospirillales","Rickettsiales","SAR11 clade","SAR202 clade", "SAR324 clade", "SAR406 clade", "Synechococcales","Other")

`%notin%` <- Negate(`%in%`)

seqdata1 <- seqdata1 %>%
  as_tibble() %>%
  arrange(Family) %>%
  relocate(Family) %>%
  mutate(Family = as.character(Family),
         Family = case_when(Family %notin% lin_focus_seq ~ "Other",
                            TRUE ~ Family),
         Family = factor(Family, levels = rev(lin_focus_seq)))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours
seq_colour <-c("grey70","green3","mediumpurple3","yellow" ,"maroon", "dodgerblue3", "tan","hotpink","lightpink", "darkseagreen1","coral", "deepskyblue", "forestgreen", "black", "mediumvioletred","gold","purple","white","red","cyan","brown4", "navy")
seqdata1.aggregated <- seqdata1 %>% group_by(Depth2, Profile, Family) %>% summarise(Abundance = sum(Abundance))
seq1_bar <- ggplot(data = seqdata1.aggregated) +
  geom_bar(aes(x = Depth2, y = Abundance, fill = Family), stat = "Identity", position = "fill") +
  ylab("ASV (Family) Fraction") + xlab("") + ggtitle("C)") +
  ## remove the grid background
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 18,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
  theme(title = element_text(color = "black", family="Arial",size = 20)) +
  scale_fill_manual(values =seq_colour) + facet_grid(~Profile) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(strip.text.x = element_text(color = "black", size = 20, face = "bold"))
seq1_bar


DH.expt.ps <- readRDS("DH.expt.rds")

## ASV should be in more than four samples and its abundance should be greater than 0.0015 of the total
DH.expt.qc = filter_taxa(DH.expt.ps, function(x) sum(x > 4) > (0.0015*length(x)), TRUE)
DH.expt = prune_samples(sample_sums(DH.expt.qc)>=1000, DH.expt.qc)

dh.expt <- subset_samples(DH.expt, Type1 == 'Expt') ## subset to the samples involved in the experiment
dh.expt <- prune_samples(sample_names(dh.expt) != 'DH36_E4', dh.expt)
## Plot bar charts for Expt at Family level for all ASVs
profile_bar <- names(sort(taxa_sums(dh.expt), decreasing=TRUE))
profile_data <- transform_sample_counts(dh.expt, function(OTU) OTU/sum(OTU))
profile_plot <- prune_taxa(profile_bar, profile_data) 
bar_plot <- tax_glom(profile_plot, "Family")
plot_bar(bar_plot, x="Sample_ID", fill ="Family") + geom_bar(stat="identity")

## Export the bar plot to ggplot to allow formating
gg <- plot_bar(bar_plot, x="Sample_ID", fill ="Family")
seqdata2 <- gg$data

#Order the factors for Day
seqdata2 <- seqdata2 %>%
  arrange(Day2) %>%
  mutate(Day2 = factor(Day2, levels = c("Day 0","Day 2","Day 6","Day 9" ,"Day 12","Day 21")))

#Order the factors for Treatment
seqdata2 <- seqdata2 %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))



`%notin%` <- Negate(`%in%`)
 ## reduce family to the more abundant ASVs using lin_focus_seq file 
seqdata2 <- seqdata2 %>%
    as_tibble() %>%
    arrange(Family) %>%
    relocate(Family) %>%
    mutate(Family = as.character(Family),
           Family = case_when(Family %notin% lin_focus_seq ~ "Other",
                              TRUE ~ Family),
           Family = factor(Family, levels = rev(lin_focus_seq)))
  
  ## To display your custom color pallets remove  +
  ## scale_y_continuous(labels = scales::percent) from the ggplot command
  seqdata2.aggregated <- seqdata2 %>% group_by(Day2, Treatment, Family) %>% summarise(Abundance = sum(Abundance))
  seq2_bar <- ggplot(data = seqdata2.aggregated) + 
   geom_bar(aes(x = Day2, y = Abundance, fill = Family), stat = "Identity", position = "fill") +
    ylab("") + xlab("") + 
    ## remove the grid background
    theme_classic() +
    ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
    theme(text = element_text(size=16, color = "black", family="Arial")) +
    theme(axis.text.x = element_text(color = "black", family="Arial", size = 18,angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", family="Arial", size = 18)) +
    scale_fill_manual(values =seq_colour) + facet_grid(~Treatment) +
    guides(fill = guide_legend(ncol = 1)) +
    theme(strip.text.x = element_text(color = "black", size = 20, face = "bold"))
  seq2_bar

  ##Combining the two ASV plots into one with one legend.
  figure5c <- ggarrange(seq1_bar, seq2_bar, ncol = 2, nrow = 1,
                        common.legend = TRUE, legend = "right",
                        align = "hv", widths = c(1,3))
  figure5c
  
  ggsave("Figure5c.png", plot = figure5c, width = 14, height = 6, dpi = 300, units = "in")

  figure5 <- ggarrange(figure5a, figure5b, figure5c, ncol = 1, nrow = 3,
                       common.legend = FALSE,
                       align = "hv")
  figure5
  
  ggsave("Figure5.png", plot = figure5, width = 12, height = 18, dpi = 300, units = "in")
  