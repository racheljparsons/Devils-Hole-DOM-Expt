## Code written by Shuting Liu and modified by Rachel Parsons in 2023 with help from Simon Biggs

library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
setwd("/Users/rachelbiggs/Library/CloudStorage/Dropbox/DH DOM Manuscript/Manuscript 2023/DH36 R Plots/AA")

aa <- read.csv(file = "AA_plot2.csv", header = TRUE)

#Order the factors for Day
aa <- aa %>%
  arrange(Day2) %>%
  mutate(Day2 = factor(Day2, levels = c("Day 0", 
                                        "Day 2",
                                        "Day 6",
                                        "Day 12")))

#Order the factors for Treatment
aa <- aa %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))



## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours
aa_colour <-c("navy", "gold", "red4", "forestgreen", "purple4", "lightblue2", "lightgoldenrod1", "lightcoral", "lightgreen", "mediumpurple1", "orange", "mediumturquoise","lightpink","lightslategrey", "violetred", "black")

ggplot(data = aa) +
  geom_bar(aes(x = Day2, y = Conc, fill = Amino_Acid), stat = "Identity", position = "fill") +
  ylab("Molar Percent") + xlab("") +
  ## remove the grid background
  theme_bw() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=20, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 16, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 16)) +
  scale_fill_manual(values = aa_colour) + facet_grid(~Treatment) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 20, face = "bold"))
aa

library(gridExtra)
library("ggpubr")

## Load the AA depth profile for DH36 and make Depth a factor
aa_profile <- read.csv(file = "AA_profile.csv", header = TRUE)
aa_profile$Depth <- as.factor(aa_profile$Depth)

##Order the factors for Depth
aa_profile <- aa_profile %>%
  arrange(Depth) %>%
  mutate(Depth = factor(Depth, levels = c("24", "23", "22", "20", "15", "10", "5", "1")))

## To display your custom color pallets remove  +
## scale_y_continuous(labels = scales::percent) from the ggplot command

## Set colours
aa_colour <-c("navy", "gold", "red4", "forestgreen", "purple4", "lightblue2", "lightgoldenrod1", "lightcoral", "lightgreen", "mediumpurple1", "orange", "mediumturquoise","lightpink","lightslategrey", "violetred", "black")

tdaa <- ggplot(data = aa_profile) +
  geom_bar(aes(x = Depth, y = Conc, fill = Amino_Acid), stat = "Identity", position = "fill") +
  ylab("Percent Molar Concentration") + xlab("Depth (m)") +
  ## remove the grid background
  theme_bw() +
  theme_classic() +
  ## increase font size, pick font, change font colour, change axis text angle and height and width adjustments
  theme(text = element_text(size=16, color = "black", family="Arial")) +
  theme(axis.text.x = element_text(color = "black", family="Arial", size = 14, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 14)) +
  scale_fill_manual(values = aa_colour) + coord_flip() + facet_grid(~Month) +
  theme(strip.text.x = element_text(color = "black", family="Arial",size = 18, face = "bold"))
tdaa

ggsave("FigureS2.png", plot = tdaa, width = 10, height = 6, dpi = 300, units = "in")
