library(tidyverse)
library(ggplot2)
library(ggtext)
library(gridExtra)

library("ggpubr")  ### use this for arranging multiple ggplots over multiple pages 
## library("here")  ## to allow you to designate file directory when saving images... see end of script.

#setwd("/Users/rachelbiggs/Library/CloudStorage/Dropbox/DH DOM Manuscript/Manuscript 1623/DH36 R Plots/Line Graphs")

line <- read.csv(file = "Line2.csv", header = TRUE)

line <- line %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

glu.graph <- ggplot(line, aes(x=Day, y=Glu))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Glu),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Glu-Glu_SD, ymax=Glu+Glu_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Glutamic Acid nmol L"^-1), colour = "Treatment") +
  ggtitle("C)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(50, 250)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
glu.graph

ser.graph <- ggplot(line, aes(x=Day, y=Ser))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Ser),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Ser-Ser_SD, ymax=Ser+Ser_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Serine nmol L"^-1), colour = "Treatment") +
  ggtitle("F)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(50, 300)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
ser.graph

gly.graph <- ggplot(line, aes(x=Day, y=Gly))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Gly),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Gly-Gly_SD, ymax=Gly+Gly_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Glycine nmol L"^-1), colour = "Treatment") +
  ggtitle("D)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(100, 800)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
gly.graph

asp.graph <- ggplot(line, aes(x=Day, y=Asp))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Asp),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Asp-Asp_SD, ymax=Asp+Asp_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Aspartic Acid nmol L"^-1), colour = "Treatment") +
  ggtitle("B)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(50, 250)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
asp.graph

ala.graph <- ggplot(line, aes(x=Day, y=Ala))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Ala),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Ala-Ala_SD, ymax=Ala+Ala_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Alanine nmol L"^-1), colour = "Treatment") +
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(50, 250)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
ala.graph

leu.graph <- ggplot(line, aes(x=Day, y=Leu))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Leu),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Leu-Leu_SD, ymax=Leu+Leu_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Leucine nmol L"^-1), colour = "Treatment") +
  ggtitle("E)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0,2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
leu.graph

tdaa.graph <- ggplot(line, aes(x=Day, y=TDAA))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$TDAA),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=TDAA-TDAA_SD, ymax=TDAA+TDAA_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("TDAA nmol L"^-1), colour = "Treatment") +
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 12)) +
  scale_y_continuous(limits = c(500, 2500)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
tdaa.graph

figureS3 <- ggarrange(ala.graph, asp.graph, glu.graph, gly.graph, leu.graph, ser.graph, ncol = 2, nrow = 3,
                      common.legend = TRUE, legend = "bottom",
                      align = "hv")
figureS3

# then save using ggsave, can specific size and resolution

ggsave("FigureS3.png", plot = figureS3, width = 8, height = 10, dpi = 300, units = "in")

