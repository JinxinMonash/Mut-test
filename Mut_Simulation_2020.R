#R version 3.6.2
"""
Created on 2020

Mutsimulation

@author: Jason

"""

#library package
library(tidyverse)
library(readxl)
library(matrixStats)
#Read the gene length data
#collect gene name and length information
#The data convert from GFF annotation file
data <- read_excel("GeneLength_genome.xlsx")

#data include gene_name, start_position, end_position, gene_description 
data <- read_excel("Zhao_simulation/GeneLength_genome0826.xlsx")
data <- data %>% 
  dplyr::select(locus_tag,cds_start,cds_end)
trial_result <- data
trial_result <- data %>% 
  mutate(Gene_Number=rep(0))

#Get the gene length and length frequency
data <- data %>% 
  mutate(length=abs(cds_end-cds_start))
data <- data %>% 
  mutate(length_freq=length/sum(data$length))

#Output file to show overall frequency distribution of gene lengths
write.csv(data,"gene_freq.csv",row.names = FALSE)

#Choosing given number of target genes by probability density function defined by gene size
#This is repeated for the selected number of simulations
#Number of simulations
trialmax <- 20000
for(i in 1:trialmax){ 
  trial <- data %>% slice_sample(weight_by = length_freq,n=99,replace=TRUE)
  trial <- as.data.frame(table(trial$locus_tag))
  colnames(trial)[1] <- "locus_tag"
  colnames(trial)[2] <- paste0("trial", i)
  trial_result <- trial_result %>% 
    left_join(trial,by="locus_tag")
  trial_result[paste0("trial", i)][is.na(trial_result[paste0("trial", i)])] <- 0
}

#Writing simulation results to CSV
write.csv(trial_result,"Simulation_Gene_20000_R_T.csv",row.names = FALSE)

#Open the simulation result
simulation_data <- read.csv("Simulation_Gene_20000.csv")
simulation_data[, 5:20004] <- sapply(simulation_data[, 5:20004], as.numeric)
#Collecting count stats
simulation_data <- simulation_data %>% 
  mutate(avg_number=rowMeans(simulation_data[,5:20004])) %>% 
  mutate(sd_number=rowSds(as.matrix(simulation_data[,5:20004]))) %>% 
  dplyr::select(locus_tag,avg_number,sd_number)
#Collect identities and information for genes of interest
interest <- read.csv("gene_interest.csv")

#Gene_name, mutation number for different condition
interest <- interest %>% 
  filter(is.na(Total)==FALSE)
interest$Control <- as.numeric(interest$Control)
interest$X0.5MIC <- as.numeric(interest$X0.5MIC)
interest$X1MIC <- as.numeric(interest$X1MIC)
interest$X4MIC <- as.numeric(interest$X4MIC)
interest$X8MIC <- as.numeric(interest$X8MIC)
interest$X16MIC <- as.numeric(interest$X16MIC)
#No mutation is 0 times
interest[is.na(interest)==TRUE] <- 0
interest <- interest %>% 
  mutate(lowcondition=X0.5MIC+X1MIC) %>% 
  mutate(highcondition=X4MIC+X8MIC+X16MIC)
interest <- interest %>% 
  left_join(simulation_data,by="locus_tag")
#Control condition gene of interest
#Low condition gene of interest
#High condition gene of interest
#Calculate the p value
interest <- interest %>% 
  mutate(TotalPvalue=
           pnorm(Total/18, mean = avg_number, sd = sd_number, lower.tail = FALSE))
interest <- interest %>% 
  mutate(ControlPvalue=
           pnorm(Control/3, mean = avg_number, sd = sd_number, lower.tail = FALSE))
interest <- interest %>% 
  mutate(LowPvalue=
           pnorm(lowcondition/6, mean = avg_number, sd = sd_number, lower.tail = FALSE))
interest <- interest %>% 
  mutate(highPvalue=
           pnorm(highcondition/9, mean = avg_number, sd = sd_number, lower.tail = FALSE))
interest <- interest %>% 
  mutate(Drug=lowcondition+highcondition)
interest <- interest %>% 
  mutate(DrugPvalue=
           pnorm(Drug/15, mean = avg_number, sd = sd_number, lower.tail = FALSE))
write.csv(interest,"interest_gene.csv",row.names = FALSE)



