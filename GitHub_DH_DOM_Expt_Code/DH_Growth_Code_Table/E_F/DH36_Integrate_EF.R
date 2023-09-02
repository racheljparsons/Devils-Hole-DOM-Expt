#####################
#DH36 Experiment
#####################

##This script was adapted to wrangle, process, and analyze the experimental data for growth in E,F treatments removing the initial death phase

##Rachel Parsons March 4th, 2023; adapted from Shuting Liu, Dec 12, 2018, adapted from Nicholas Huynh, November 6, 2018

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
setwd("/Users/rachelbiggs/Library/CloudStorage/Dropbox/DH DOM Manuscript/R code/Integration Code/DH_Growth_Code_Table/E_F")
ef.df <- read.csv("DH36_Growth_EF.csv")

#Part 3. Add derived variables ----

#add a column of the natural logarithm of the cell counts

ef.df <- 
ef.df %>% 
  mutate(lncellsL = round(log(Cells_L), 2))

#Part 3a. Interpolate the organic carbon values ----

ef.df <- 
  ef.df %>%
  group_by(Expt, Sample) %>%
  mutate(INT_DOC = na.approx(DOC,Day, na.rm=F,rule=2)) %>%
ungroup()
  
#######GROWTH CURVER FOR BACT AbUNDANCE#####

#Part 4. Fitting the growth curves to a logistic equation ----

#Here, we are using the growthcurver package to fit growth curve data to the standard form of the logistic equation and then return a data table with population-level information. Specifically, the carrying capacity k, specific growth rate r, the initial population size N0, the goodness of fit to the logistic equation sigma as well as the degrees of freedom df, the time at which 1/2 carrying capcicity is reached t_mid, the doubling time t_gen, the area under the logistic curve auc_l, and the emprical area under the experimental curve auc_e. 

#4a. Prepare data for growthcurver ----

gcurve_ef.df <- 
 ef.df %>% 
  group_by(Expt, Sample) %>%
  #initial cells_ml value for each experiment will be subtracted from the cells_ml values 
  mutate(delta_cells = Cells_L - Cells_L[which.min(Day)],
         delta_doc = INT_DOC[which.min(Day)] - INT_DOC,) %>% 
  ungroup() 

write.csv(gcurve_ef.df,"gcurve_DH_EF.csv")

#Exclude the death phase in E and F Treatment D/S
gcurve_ef_nodeath.df <- gcurve_ef.df [-c(1:8,23:24),]  ##Remove the death phase of the growth curve
#split the dataframe above by sample
gcurve_ef_nodeath.list <- split(gcurve_ef_nodeath.df, gcurve_ef.df$Sample) 
#store the names of each list object (i.e. ID2)
headers = names(gcurve_ef_nodeath.list)

#4b. Create a function/for loop that plots each of the curves ----

#first apply the summarize growth function to each of the experiments
gcurveplot.func <- function(elf){
  gc_ef_fit <- SummarizeGrowth(elf$Day, elf$delta_cells) #can use t_trim to trim out death phase if did not remove before
}
gcplot_ef.list <- lapply(gcurve_ef_nodeath.list, gcurveplot.func)

#4c. save the plots as a pdf ----
pdf("growthcurves_DH_EF.pdf")
for (i in 1:length(gcplot_ef.list)) {
  plot(gcplot_ef.list[[i]], main = names(gcplot_ef.list[i]) ) 
}
dev.off()

#4d. Create a function that returns growthcurver output fromt the test curves as a dataframe ----

gcurve.func <- function(narwhal){
  gc.fit <- SummarizeGrowth(narwhal$Day, narwhal$delta_cells)
  gc_output_ef.df  <- as.data.frame(unlist(gc.fit$vals), stringsAsFactors = FALSE)
}

#4e. Apply the function to all the growth curves ----
gc_out_ef.list <- lapply(gcurve_ef_nodeath.list, gcurve.func)
#save the list as a data frame
gc_out_ef.df <- data.frame(gc_out_ef.list, stringsAsFactors = FALSE)
#transpose the data frame
gc_out_ef.df <- as.data.frame(t(gc_out_ef.df), stringsAsFactors = FALSE)
#replace the character strings
#gc_out.df$note <- gsub("cannot fit data", "1", gc_out.df$note)
#gc_out.df$note <- gsub("questionable fit", "2", gc_out.df$note)
#coerce the data frame to be one of numerics, not characters
gc_out_ef.df <- as.data.frame(sapply(gc_out_ef.df, as.numeric)) 
#reapply the sample names as the row names
rownames(gc_out_ef.df) <- headers
#make the row names the first column of the data frame and change the column name to "ID2"
gc_out_ef.df <-setDT(gc_out_ef.df, keep.rownames = TRUE)[] #make the rownames into the first column
colnames(gc_out_ef.df)[1] <- "Bottle"

#4f. Tidy the growth curve data and find the transition points of the growth phases----

gc_out_ef.df <- 
  gc_out_ef.df %>%
  mutate_at(vars(k, k_se, n0, n0_se, sigma, auc_l, auc_e),funs(round(., 0))) %>% #auc_l area under curve from fitting, auc_e from measurement, these are auc from 0 to end of time, k carrying capacity,se standard error,p value,n0 initial population size, r growth rate,tmid time to reach half carrying capacity, sigma residual standard error from fit,t_gen doubling time
  mutate_at(vars(t_mid, t_gen),funs(round(., 1))) %>%
  mutate_at(vars(k_p, n0_p, r, r_p),funs(round(., 2))) %>% 
  mutate(stationary = t_mid*2)  #stationary phase defined as 2 times half the time it takes reach the carrying capacity (delta).
  
write.csv(gc_out_ef.df,"gc_out_DH_EF.csv")

gc_out_select_ef.df <- gc_out_ef.df %>%
  select(Bottle, stationary,t_gen)


#4g. merge the growthcurve output with the input dataframe

gcurve_ef.df <- 
  gcurve_ef.df %>% 
  cross_join(., gc_out_select_ef.df)

write.csv(gcurve_ef.df,"gcurve_new_DH_EF.csv")

#there we will modify the data frame in excel, adding the stationary time points as new rows and then re-import the file. A,D,J pretty close to real time points, don't need interpolation
#we will then interpolate data for those missing points
#this is necessary to be able to calculate accurate areas under the curves and thus, bge
########IMPORT DATA##########

int.df <- read.csv("gcurve_check_DH_EF.csv")

int.df <- int.df %>% 
  select(Expt, Timepoint, Date, Sample, Treatment,Day,stationary,t_gen, Cells_L, DOC, lncellsL, INT_DOC, delta_cells, delta_doc) %>% 
  group_by(Expt, Sample) %>%
  mutate(INTERP_cells = na.approx(Cells_L, Day, na.rm=F),
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

write.csv(auc_int.df,"gauc_int_DH_EF.csv")
auc_int.df <- read.csv("gauc_int_DH_EF.csv")

#add a new column to each data frame with the specific growth rate
#growth rates are calculated as the slope of the exponential phase for each bottle

e.df <- filter(auc_int.df, Sample == "E")
f.df <- filter(auc_int.df, Sample == "F")

e.df <- mutate(e.df, µ = round(lm(lncellsL~Day,e.df[c(3:7),])$coeff[[2]],2))
f.df <- mutate(f.df, µ = round(lm(lncellsL~Day,f.df[c(2:7),])$coeff[[2]],2))

#put all of the dataframes back together
DHgrowth <- rbind(e.df, f.df)

DH_calc <- DHgrowth %>%
  select(Expt, Sample, Day, stationary, t_gen, DOC,INT_DOC, delta_cells, delta_doc, INTERP_cells, INTERP_DOC, INTERP_delta_cells,INTERP_delta_doc,AUC_cells, AUC_doc,AUC_cells_norm,AUC_doc_norm, µ) 


#save the data frame as a csv file 
write.csv(DH_calc, "DH_Calculated_Master_EF.csv")

