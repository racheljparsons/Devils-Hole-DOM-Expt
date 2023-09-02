#####################
#DH36 Experiment
#####################

#This script was written to wrangle, process, and analyze the experimental data (BC using boivolume)

#Rachel Parsons March 4th, 2023; adapted from Shuting Liu, Dec 12, 2018, adapted from Nicholas Huynh, November 6, 2018

######PREP DATA########

#Part 1. load necessary packages and data ----
library(tidyverse)
library(data.table)
library(growthcurver)
library(oce)
library(zoo)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ggplot2)

#Part 2. Reading in the dataset ----
setwd("/Users/rachelbiggs/Library/CloudStorage/Dropbox/DH DOM Manuscript/R code/Integration Code/DH_Growth_Code_Table")
master.df <- read.csv("DH36_Growth.csv")

#Part 3. Add derived variables ----

#add a column of the natural logarithm of the cell counts

master.df <- 
  master.df %>% 
  mutate(lncellsml = round(log(Cells_mL), 2))

#Part 3a. Interpolate the organic carbon values ----

master.df <- 
  master.df %>%
  group_by(Expt, Sample) %>%
  mutate(INT_DOC = na.approx(DOC,Day, na.rm=F,rule=2)) %>%
ungroup()
  
#######GROWTH CURVER FOR BACT AbUNDANCE#####

#Part 4. Fitting the growth curves to a logistic equation ----

#Here, we are using the growthcurver package to fit growth curve data to the standard form of the logistic equation and then return a data table with population-level information. Specifically, the carrying capacity k, specific growth rate r, the initial population size N0, the goodness of fit to the logistic equation sigma as well as the degrees of freedom df, the time at which 1/2 carrying capcicity is reached t_mid, the doubling time t_gen, the area under the logistic curve auc_l, and the emprical area under the experimental curve auc_e. 

#4a. Prepare data for growthcurver ----

gcurve.df <- 
 master.df %>% 
  group_by(Expt, Sample) %>%
  #initial cells_ml value for each experiment will be subtracted from the cells_ml values 
  mutate(delta_cells = Cells_mL - Cells_mL[which.min(Day)],
         delta_doc = INT_DOC[which.min(Day)] - INT_DOC,) %>% 
  ungroup() 

write.csv(gcurve.df,"gcurve_DH.csv")

#Exclude the death phase in E and F Treatment D/S
gcurve_nodeath.df <- gcurve.df [-c(5:6,11:12,59:60,65:66,69:84),]  ##From Shuting: this is where I exclude the different death phase from your original code


#split the dataframe above by sample
gcurve_nodeath.list <- split(gcurve_nodeath.df, gcurve_nodeath.df$Sample) ###From Shuting: this should be gcurve_nodeath not gcurve.df
#store the names of each list object (i.e. ID2)
headers = names(gcurve_nodeath.list)

#4b. Create a function/for loop that plots each of the curves ----

#first apply the summarize growth function to each of the experiments
gcurveplot.func <- function(elf){
  gc_fit <- SummarizeGrowth(elf$Day, elf$delta_cells) #can use t_trim to trim out death phase if did not remove before
}
gcplot.list <- lapply(gcurve_nodeath.list, gcurveplot.func)

#4c. save the plots as a pdf ----
pdf("growthcurves_DH.pdf")
for (i in 1:length(gcplot.list)) {
  plot(gcplot.list[[i]], main = names(gcplot.list[i]) ) 
}
dev.off()

#4d. Create a function that returns growthcurver output fromt the test curves as a dataframe ----

gcurve.func <- function(narwhal){
  gc.fit <- SummarizeGrowth(narwhal$Day, narwhal$delta_cells)
  gc_output.df  <- as.data.frame(unlist(gc.fit$vals), stringsAsFactors = FALSE)
}

#4e. Apply the function to all the growth curves ----
gc_out.list <- lapply(gcurve_nodeath.list, gcurve.func)
#save the list as a data frame
gc_out.df <- data.frame(gc_out.list, stringsAsFactors = FALSE)
#transpose the data frame
gc_out.df <- as.data.frame(t(gc_out.df), stringsAsFactors = FALSE)
#replace the character strings
#gc_out.df$note <- gsub("cannot fit data", "1", gc_out.df$note)
#gc_out.df$note <- gsub("questionable fit", "2", gc_out.df$note)
#coerce the data frame to be one of numerics, not characters
gc_out.df <- as.data.frame(sapply(gc_out.df, as.numeric)) 
#reapply the sample names as the row names
rownames(gc_out.df) <- headers
#make the row names the first column of the data frame and change the column name to "ID2"
gc_out.df <-setDT(gc_out.df, keep.rownames = TRUE)[] #make the rownames into the first column
colnames(gc_out.df)[1] <- "Bottle"

#4f. Tidy the growth curve data and find the transition points of the growth phases----

gc_out.df <- 
  gc_out.df %>%
  mutate_at(vars(k, k_se, n0, n0_se, sigma, auc_l, auc_e),funs(round(., 0))) %>% #auc_l area under curve from fitting, auc_e from measurement, these are auc from 0 to end of time, k carrying capacity,se standard error,p value,n0 initial population size, r growth rate,tmid time to reach half carrying capacity, sigma residual standard error from fit,t_gen doubling time
  mutate_at(vars(t_mid, t_gen),funs(round(., 1))) %>%
  mutate_at(vars(k_p, n0_p, r, r_p),funs(round(., 2))) %>% 
  mutate(stationary = t_mid*2)  #stationary phase defined as 2 times half the time it takes reach the carrying capacity (delta).
  
write.csv(gc_out.df,"gc_out_DH.csv")

gc_out_select.df <- gc_out.df %>%
  select(Bottle, stationary,t_gen)


#4g. merge the growthcurve output with the input dataframe

gcurve.df <- 
  gcurve.df %>% 
  cross_join(., gc_out_select.df)

write.csv(gcurve.df,"gcurve_new_DH.csv")

#there we will modify the data frame in excel, adding the stationary time points as new rows and then re-import the file. A,D,J pretty close to real time points, don't need interpolation
#we will then interpolate data for those missing points
#this is necessary to be able to calculate accurate areas under the curves and thus, bge
########IMPORT DATA##########

int.df <- read.csv("gcurve_check_DH.csv")

int.df <- int.df %>% 
  select(Expt, Timepoint, Date, Sample, Treatment,Day,stationary,t_gen, Cells_mL, DOC, lncellsml, INT_DOC, delta_cells, delta_doc) %>% 
  group_by(Expt, Sample) %>%
  mutate(INTERP_cells = na.approx(Cells_mL, Day, na.rm=F),
         INTERP_DOC = na.approx(INT_DOC, Day, na.rm=F,rule=2),
         INTERP_delta_cells=na.approx(delta_cells,Day, na.rm=F),
         INTERP_delta_doc = na.approx(delta_doc, Day,na.rm = F))

#Part 5. calc int area under the curve ---- using true measured data
auc_int.df <- 
  int.df %>% 
  group_by(Expt, Sample) %>%
  mutate(
    AUC_cells = integrateTrapezoid(Day, INTERP_delta_cells, type="cA"),
    AUC_doc = integrateTrapezoid(Day, INTERP_delta_doc, type="cA"), #delta_doc already include interpolated all time points data
  ) %>% 
  distinct()     #unique rows  
#this gives Integrated values for cells and DOC for every time points

#add column for AUC normalized to day
auc_int.df <- 
  auc_int.df %>% 
  group_by(Expt, Sample) %>%
  mutate(
    AUC_cells_norm = AUC_cells/Day,
    AUC_doc_norm = AUC_doc/Day) 

write.csv(auc_int.df,"gauc_int_DH.csv")
auc_int.df <- read.csv("gauc_int_DH.csv")

####PLOTTING####
gcurve.df <- read.csv("gcurve_DH.csv")
doc.graph <- ggplot(gcurve.df, aes(x=Day, y=DOC)) +
  geom_point(aes(colour=Sample, size =1)) +
  geom_line(data=gcurve.df[!is.na(gcurve.df$DOC),],aes(colour=Sample, size=0.8)) +
  guides(size=FALSE) +
  labs(x="Days", y="DOC (µM C)", colour = "Sample") +
  ggtitle("DOC") +
  theme_classic() +
  scale_color_manual(values=c('dodgerblue3','dodgerblue','seagreen','seagreen2', 'orangered','orange')) +
  scale_x_continuous(breaks = seq(min(gcurve.df$Day), max(gcurve.df$Day), by = 10)) + 
  ## scale_y_continuous(breaks = seq(min(gcurve.df$DOC), max(gcurve.df$DOC), by = 1)) + 
  theme(title = element_text(size =17),
        axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.title=element_text(size=16),  
        legend.text = element_text(size=16))
doc.graph

## plot ln of cell growth need to add colour and error
lncells.graph <- ggplot(gcurve.df, aes(x=Day, y=lncellsml)) +
  geom_point(aes(colour=Sample, size =1)) +
  geom_line(aes(colour=Sample, size=0.8)) +
  guides(size=FALSE) +
  labs(x="Days", y="ln(Cells/ml)", colour = "Sample") +
  ggtitle("ln(Cells/ml)") +
  theme_classic() +
  scale_color_manual(values=c('dodgerblue3','dodgerblue','seagreen','seagreen2', 'orangered','orange')) +
  scale_x_continuous(breaks = seq(min(gcurve.df$Day), max(gcurve.df$Day), by = 10)) + 
  scale_y_continuous(breaks = seq(min(gcurve.df$lncellsml), max(gcurve.df$lncellsml), by = 0.2)) + 
  theme(title = element_text(size =17),
        axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.title=element_text(size=16),  
        legend.text = element_text(size=16))
lncells.graph

#add a new column to each data frame with the specific growth rate
#growth rates are calculated as the slope of the exponential phase for each bottle

a.df <- filter(auc_int.df, Sample == "A") 
b.df <- filter(auc_int.df, Sample == "B") 
c.df <- filter(auc_int.df, Sample == "C")
d.df <- filter(auc_int.df, Sample == "D")
e.df <- filter(auc_int.df, Sample == "E")
f.df <- filter(auc_int.df, Sample == "F")


a.df <- mutate(a.df, µ = round(lm(lncellsml~Day,a.df[c(5:12),])$coeff[[2]],2))
b.df <- mutate(b.df, µ = round(lm(lncellsml~Day,b.df[c(4:12),])$coeff[[2]],2))
c.df <- mutate(c.df, µ = round(lm(lncellsml~Day,c.df[c(3:9),])$coeff[[2]],2))
d.df <- mutate(d.df, µ = round(lm(lncellsml~Day,d.df[c(3:9),])$coeff[[2]],2))
e.df <- mutate(e.df, µ = round(lm(lncellsml~Day,e.df[c(5:9),])$coeff[[2]],2))
f.df <- mutate(f.df, µ = round(lm(lncellsml~Day,f.df[c(5:9),])$coeff[[2]],2))

#put all of the dataframes back together
DHgrowth <- rbind(a.df, b.df, c.df, d.df, e.df, f.df)

DH_calc <- DHgrowth %>%
  select(Expt, Sample, Day, stationary, t_gen, DOC,INT_DOC, delta_cells, delta_doc, INTERP_cells, INTERP_DOC, INTERP_delta_cells,INTERP_delta_doc,AUC_cells, AUC_doc,AUC_cells_norm,AUC_doc_norm, µ) 


#save the data frame as a csv file 
write.csv(DH_calc, "DH_Calculated_Master2.csv")
         