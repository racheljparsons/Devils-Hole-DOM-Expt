### Tukey code written by Rachel Parsons in 2023

# Load required libraries
library(tidyr)
library(stats)

tukey <- read.csv(file = "DH_Tukey.csv", header = TRUE)

# Perform Tukey test 
tukey_bact <- TukeyHSD(aov(Bacteria ~ Treatment, data = tukey))
# Print the Tukey test results
print(tukey_bact)

tukey_doc <- TukeyHSD(aov(DOC ~ Treatment, data = tukey))
# Print the Tukey test results
print(tukey_doc)

tukey_gene <- read.csv(file = "Gene.csv", header = TRUE)

tukey_cbbl <- TukeyHSD(aov(cbbL ~ Treatment, data = tukey))
# Print the Tukey test results
print(tukey_cbbl)

tukey_amoA <- TukeyHSD(aov(amoA ~ Treatment, data = tukey_gene))
# Print the Tukey test results
print(tukey_amoA)

tukey_line <- read.csv(file = "Line2.csv", header = TRUE)

tukey_TDAA <- TukeyHSD(aov(TDAA_C ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_TDAA)

tukey_yield <- TukeyHSD(aov(Yield ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_yield)

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

tukey_asp <- TukeyHSD(aov(Asp ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_asp)

tukey_ala <- TukeyHSD(aov(Ala ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_ala)

tukey_glu <- TukeyHSD(aov(Glu ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_glu)

tukey_ser <- TukeyHSD(aov(Ser ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_ser)

tukey_leu <- TukeyHSD(aov(Leu ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_leu)

tukey_gly <- TukeyHSD(aov(Gly ~ Treatment, data = tukey_line))
# Print the Tukey test results
print(tukey_gly)

diff <- read.csv(file = "Diff2.csv", header = TRUE)
diff2 <- diff[diff$Day == '6 Days',]
diff3 <- diff[diff$Day == '2 Days',]

# Perform Tukey test on Bacteria difference data
diff2_bact <- TukeyHSD(aov(Bacteria ~ Treatment, data = diff2))
# Print the Tukey test results
print(diff2_bact)

# Perform Tukey test 
diff3_doc <- TukeyHSD(aov(DOC ~ Treatment, data = diff3))
# Print the Tukey test results
print(diff3_doc)

tukey2 <- tukey[tukey$Timepoint == "T7",]
tukey3 <- tukey[tukey$Timepoint == "T8",]
tukey4 <- tukey[tukey$Timepoint == "T4",]
tukey5 <- tukey[tukey$Timepoint == "T13",]
tukey6 <- tukey[tukey$Timepoint == "T11",]

tukey2_bact <- TukeyHSD(aov(Bacteria ~ Treatment, data = tukey2))
# Print the Tukey test results
print(tukey2_bact)

tukey3_bact <- TukeyHSD(aov(Bacteria ~ Treatment, data = tukey3))
# Print the Tukey test results
print(tukey3_bact)

tukey3_di <- TukeyHSD(aov(DI ~ Treatment, data = tukey3))
# Print the Tukey test results
print(tukey3_di)

tukey4_cbbL <- TukeyHSD(aov(cbbL ~ Treatment, data = tukey4))
# Print the Tukey test results
print(tukey4_cbbL)

tukey5_sar202 <- TukeyHSD(aov(SAR202 ~ Treatment, data = tukey5))
# Print the Tukey test results
print(tukey5_sar202)

tukey4_syn <- TukeyHSD(aov(Synecococcus ~ Treatment, data = tukey4))
# Print the Tukey test results
print(tukey4_syn)

tukey6_nh4 <- TukeyHSD(aov(Ammonia ~ Treatment, data = tukey6))
# Print the Tukey test results
print(tukey6_nh4)

tukey5_no2 <- TukeyHSD(aov(Nitrite ~ Treatment, data = tukey5))
# Print the Tukey test results
print(tukey5_no2)

tukey5_amoA <- TukeyHSD(aov(amoA ~ Treatment, data = tukey5))
# Print the Tukey test results
print(tukey5_amoA)

tukey5_thaum <- TukeyHSD(aov(Thaumarcheota ~ Treatment, data = tukey5))
# Print the Tukey test results
print(tukey5_thaum)