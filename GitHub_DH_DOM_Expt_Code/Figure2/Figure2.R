## Code writen in 2023 by Rachel Parsons with assistance from Simon Biggs

library(tidyverse)
library(ggplot2)
library(ggtext)
library(gridExtra)

library("ggpubr")  ### use this for arranging multiple ggplots over multiple pages 
## library("here")  ## to allow you to designate file directory when saving images... see end of script.

#setwd("/Users/rachelbiggs/Library/CloudStorage/Dropbox/DH DOM Manuscript/Manuscript 1623/DH36 R Plots/Line Graphs")

setwd("~/Desktop/Rachel_R_example")
line <- read.csv(file = "Line2.csv", header = TRUE)

line <- line %>%
  arrange(Treatment) %>%
  mutate(Treatment = factor(Treatment, levels = c("S/S", 
                                                  "S/D",
                                                  "D/S")))

bact.graph <- ggplot(line, aes(x=Day, y=Bacteria))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line,aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Bacteria-Bacteria_SD, ymax=Bacteria+Bacteria_SD, colour=Treatment), width=0.5) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("Prokaryotic Abundance (x10" ^9 ~ "cells L"^-1~")"), colour = "Treatment") + # use expression() function, have to kept blocks of text in "" then ^ to start superscript and ~ to stop
  # labs(x="Days", y=expression("Prokaryotes (x10" [9] ~ "cells L"[-1] ~")"), colour = "Treatment")) ### the [] would produce a subscript
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(5, 25)) +
  theme(title = element_text(size=16, family="Times New Roman"),
        axis.title = element_text(size=16, family="Times New Roman"),
        axis.text = element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16, family="Times New Roman"),  
        legend.text = element_text(size=16, family="Times New Roman"))
bact.graph

doc.graph <- ggplot(line, aes(x=Day, y=DOC))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$DOC),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=DOC-DOC_SD, ymax=DOC+DOC_SD, colour=Treatment), width=0.5) +
  guides(size=FALSE) +
  labs(x="Days", y=expression("DOC µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("B)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(80, 110)) +
  theme(title = element_text(size=16, family="Times New Roman"),
        axis.title = element_text(size=16, family="Times New Roman"),
        axis.text = element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16, family="Times New Roman"),  
        legend.text = element_text(size=16, family="Times New Roman"))
doc.graph

ammonium.graph <- ggplot(line, aes(x=Day, y=Ammonium))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Ammonium),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Ammonium-NH4_SD, ymax=Ammonium+NH4_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Ammonium µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("C)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  theme(title = element_text(size=16, family="Times New Roman"),
        axis.title = element_text(size=16, family="Times New Roman"),
        axis.text = element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16, family="Times New Roman"),  
        legend.text = element_text(size=16, family="Times New Roman"))
ammonium.graph

nitrite.graph <- ggplot(line, aes(x=Day, y=Nitrite))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Nitrite),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Nitrite-NO2_SD, ymax=Nitrite+NO2_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Nitrite µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("D)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 5)) +
  theme(title = element_text(size=16, family="Times New Roman"),
        axis.title = element_text(size=16, family="Times New Roman"),
        axis.text = element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16, family="Times New Roman"),  
        legend.text = element_text(size=16, family="Times New Roman"))
nitrite.graph


## can summarise plots using common legend using ggarrange in ggpubr package

figure2 <- ggarrange(bact.graph, doc.graph, ammonium.graph, nitrite.graph, ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "bottom",
                     align = "hv")
figure2

# then save using ggsave, can specific size and resolution


ggsave("Figure2.png", plot = figure2, width = 8, height = 6, dpi = 300, units = "in")

