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

tdaa_c.graph <- ggplot(line, aes(x=Day, y=TDAA_C))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$TDAA_C),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=TDAA_C-TDAA_C_SD, ymax=TDAA_C+TDAA_C_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("TDAA C nmol L"^-1), colour = "Treatment") +
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(1000, 9000),breaks = c(1000,3000,5000,7000, 9000)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
tdaa_c.graph

yield.graph <- ggplot(line, aes(x=Day, y=Yield))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Yield),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Yield-Yield_SD, ymax=Yield+Yield_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y="% TDAA Yield", colour = "Treatment") +
  ggtitle("C)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(2, 9), breaks = c(2,3,4,5,6,7,8,9)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
yield.graph

di.graph <- ggplot(line, aes(x=Day, y=DI))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$DI),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=DI-DI_SD, ymax=DI+DI_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y="Degradation Index", colour = "Treatment") +
  ggtitle("B)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0,1,2,3)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
di.graph


gene <- read.csv(file = "Gene.csv", header = TRUE)

gene <- gene %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

library(scales)

cbbL.graph <- ggplot(gene, aes(x=Day,y=cbbl))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=gene[!is.na(gene$cbbl),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=cbbl-cbbl_SD, ymax=cbbl+cbbl_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("cbbL (x10"^3 ~"gcn/ng DNA)"), colour = "Treatment") +
  ggtitle("D)") +
  theme_classic() +
  scale_colour_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_continuous(limits = c(0, 21)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title.y = element_text(size=16, family="Arial", vjust = 0),
        axis.text.x = element_text(color = "black", family="Arial", size = 14, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", family="Arial", size = 14),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
cbbL.graph

lineage <- read.csv(file = "Lineage.csv", header = TRUE)

lineage <- lineage %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

sar202.graph <- ggplot(lineage, aes(x=Day, y=SAR202))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=lineage[!is.na(lineage$SAR202),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=SAR202-SAR202_SD, ymax=SAR202+SAR202_SD, colour=Treatment), width=0.5) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("SAR202 (x10" ^7 ~ "cells L"^-1~")"), colour = "Treatment") + # use expression() function, have to kept blocks of text in "" then ^ to start superscript and ~ to stop
  # labs(x="Days", y=expression("Prokaryotes (x10" [9] ~ "cells L"[-1] ~")"), colour = "Treatment")) ### the [] would produce a subscript
  ggtitle("E)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
sar202.graph

chlorobi.graph <- ggplot(lineage, aes(x=Day, y=Chlorobi))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=lineage[!is.na(lineage$Chlorobi),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Chlorobi-Chlorobi_SD, ymax=Chlorobi+Chlorobi_SD, colour=Treatment), width=0.5) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("Chlorobiaceae (x10" ^7 ~ "cells L"^-1~")"), colour = "Treatment") + # use expression() function, have to kept blocks of text in "" then ^ to start superscript and ~ to stop
  # labs(x="Days", y=expression("Prokaryotes (x10" [9] ~ "cells L"[-1] ~")"), colour = "Treatment")) ### the [] would produce a subscript
  ggtitle("F)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
chlorobi.graph

ggsave("Chlorobi.png", plot = chlorobi.graph, width = 6, height = 4, dpi = 300, units = "in")


syn.graph <- ggplot(lineage, aes(x=Day, y=Syn))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=lineage[!is.na(lineage$Syn),],aes(colour=Treatment), linewidth=1) +
  geom_errorbar(aes(ymin=Syn-Syn_SD, ymax=Syn+Syn_SD, colour=Treatment), width=0.5) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("Synechococcus (x10" ^7 ~ "cells L"^-1~")"), colour = "Treatment") + # use expression() function, have to kept blocks of text in "" then ^ to start superscript and ~ to stop
  # labs(x="Days", y=expression("Prokaryotes (x10" [9] ~ "cells L"[-1] ~")"), colour = "Treatment")) ### the [] would produce a subscript
  ggtitle("F)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 15), breaks = c(2.5,5,7.5,10,12.5,15)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=14, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=14, family="Arial"))
syn.graph

## can summarise plots using common legend using ggarrange in ggpubr package

figure3 <- ggarrange(tdaa_c.graph, di.graph, yield.graph, cbbL.graph, sar202.graph, syn.graph, ncol = 2, nrow = 3,
                      common.legend = TRUE, legend = "bottom",
                      align = "hv")
figure3

# then save using ggsave, can specific size and resolution

ggsave("Figure3new.png", plot = figure3, width = 8, height = 10, dpi = 300, units = "in")

