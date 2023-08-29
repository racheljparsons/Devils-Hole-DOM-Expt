# Load required libraries
library(tidyr)
library(stats)

tukey <- read.csv(file = "Tukey.csv", header = TRUE)

# Perform Tukey test
tukey_bact <- TukeyHSD(aov(Bacteria ~ Treatment, data = tukey))
# Print the Tukey test results
print(tukey_bact)

tukey_doc <- TukeyHSD(aov(DOC ~ Treatment, data = tukey))
# Print the Tukey test results
print(tukey_doc)

tukey_gene <- read.csv(file = "Gene.csv", header = TRUE)

tukey_cbbl <- TukeyHSD(aov(cbbl ~ Treatment, data = tukey_gene))
# Print the Tukey test results
print(tukey_cbbl)

tukey_amoA <- TukeyHSD(aov(amoA ~ Treatment, data = tukey_gene))
# Print the Tukey test results
print(tukey_amoA)

tukey_line <- read.csv(file = "Line2.csv", header = TRUE)

tukey_DI <- TukeyHSD(aov(DI ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_DI)

tukey_ammonia <- TukeyHSD(aov(Ammonium ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_ammonia)

tukey_nitrite <- TukeyHSD(aov(Nitrite ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_nitrite)

tukey_nitrate <- TukeyHSD(aov(Nitrate ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_nitrate)

tukey_silicate <- TukeyHSD(aov(Silicate ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_silicate)

tukey_phosphate <- TukeyHSD(aov(Phosphate ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_phosphate)

tukey_lineage <- read.csv(file = "Lineage.csv", header = TRUE)

tukey_thaum <- TukeyHSD(aov(Thaum ~ Treatment, data = tukey_lineage))
# Print the Tukey test results
print(tukey_thaum)

tukey_sar202 <- TukeyHSD(aov(SAR202 ~ Treatment, data = tukey_lineage))
# Print the Tukey test results
print(tukey_sar202)

tukey_syn <- TukeyHSD(aov(Syn ~ Treatment, data = tukey_lineage))
# Print the Tukey test results
print(tukey_syn)

tukey_chlorobi <- TukeyHSD(aov(Chlorobi ~ Treatment, data = tukey_lineage))
# Print the Tukey test results
print(tukey_chlorobi)