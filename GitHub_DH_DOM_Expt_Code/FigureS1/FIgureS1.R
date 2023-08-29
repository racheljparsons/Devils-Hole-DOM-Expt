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

toc.graph <- ggplot(line, aes(x=Day, y=TOC))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$TOC),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=TOC-TOC_SD, ymax=TOC+TOC_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("TOC µmol L"^-1), colour = "Treatment") +
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(80, 110)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
toc.graph


phosphate.graph <- ggplot(line, aes(x=Day, y=Phosphate))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Phosphate),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Phosphate-PO4_SD, ymax=Phosphate+PO4_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Phosphate µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("B)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
phosphate.graph

silicate.graph <- ggplot(line, aes(x=Day, y=Silicate))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Silicate),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Silicate-SO3_SD, ymax=Silicate+SO3_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Silicate µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("C)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 5)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
silicate.graph

nitrate.graph <- ggplot(line, aes(x=Day, y=Nitrate))  +
  geom_point(aes(colour=Treatment), size=3) +
  geom_line(data=line[!is.na(line$Nitrate),],aes(colour=Treatment), size=1) +
  geom_errorbar(aes(ymin=Nitrate-NO3_SD, ymax=Nitrate+NO3_SD, colour=Treatment), width=0.5) +
  labs(x="Days", y=expression("Nitrate µmol L"^-1), colour = "Treatment") + # use expression() then text in "", ^ to start superscript then the brackets close superscript function as no more text in Y label.
  ggtitle("A)") +
  theme_classic() +
  scale_color_manual(values=c("dodgerblue3", "mediumseagreen", "orangered")) +
  scale_x_continuous(limits = c(0, 21)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme(title = element_text(size=16, family="Arial"),
        axis.title = element_text(size=16, family="Arial"),
        axis.text = element_text(size=16, family="Arial"),
        legend.title=element_text(size=16, family="Arial"),  
        legend.text = element_text(size=16, family="Arial"))
nitrate.graph

## can summarise plots using common legend using ggarrange in ggpubr package

figureS1 <- ggarrange(nitrate.graph, phosphate.graph, silicate.graph, ncol = 2, nrow = 2,
                      common.legend = TRUE, legend = "bottom",
                      align = "hv")
figureS1

# then save using ggsave, can specific size and resolution

ggsave("FigureS1.png", plot = figureS1, width = 8, height = 8, dpi = 300, units = "in")

