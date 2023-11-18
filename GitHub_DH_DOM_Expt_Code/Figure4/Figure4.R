library(tidyverse)
library(ggplot2)
library(ggtext)
library(gridExtra)
library("ggpubr")

library(scales)

line <- read.csv(file = "Line2.csv", header = TRUE)

line <- line %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))
ammonium.graph <- ggplot(line, aes(x=Day, y=Ammonium))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Ammonium),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Ammonium-NH4_SD, ymax=Ammonium+NH4_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Total Ammonium µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
ammonium.graph

nitrite.graph <- ggplot(line, aes(x=Day, y=Nitrite))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Nitrite),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Nitrite-NO2_SD, ymax=Nitrite+NO2_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Nitrite µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("B)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 5)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
nitrite.graph

gene <- read.csv(file = "Gene.csv", header = TRUE)

gene <- gene %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))
amoA.graph <- ggplot(gene, aes(x=Day,y=amoA))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=gene[!is.na(gene$amoA),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=amoA-amoA_SD, ymax=amoA+amoA_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("amoA (x10"^3 ~"gcn/ng DNA)"), colour = "Treatment") +
  ggtitle("C)") +
  theme_classic() +
  scale_colour_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_y_continuous(limits = c(0, 15)) + 
  scale_x_continuous(limits = c(0, 21)) +
  theme(title = element_text(size=18, family="Arial"),
        axis.title.y = element_text(size=18, family="Arial", vjust = 0),
        axis.text.x = element_text(color = "black", family="Arial", size = 16, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 16),
        legend.title=element_text(size=18, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
amoA.graph

lineage <- read.csv(file = "Lineage.csv", header = TRUE)

lineage <- lineage %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

thaum.graph <- ggplot(lineage, aes(x=Day, y=Thaum))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=lineage[!is.na(lineage$Thaum),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Thaum-Thaum_SD, ymax=Thaum+Thaum_SD, colour=Treatment), width=0.2) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("Thaumarcheota (x10" ^7 ~ "cells L"^-1~")"), colour = "Treatment") +
  ggtitle("D)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 30),breaks = c(0,5,10,15,20,25,30)) + 
  theme(title = element_text(size=18, family="Arial"),
        axis.title = element_text(size=18, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=18, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
thaum.graph

figure4 <- ggarrange(ammonium.graph, nitrite.graph, amoA.graph, thaum.graph, ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "bottom",
                     align = "hv")
figure4

ggsave("Figure4.png", plot = figure4, width = 8, height = 8, dpi = 300, units = "in")