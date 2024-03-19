rm(list = ls())
library(phyloseq)
IntData = readRDS("IntegratedBacteriaData.rds")
size<-min(sample_sums(IntData))
#15033
IntData<-rarefy_even_depth(IntData,rngseed=1,sample.size=15033,replace=T)
for(Exp in c("Exp1","Exp2")){
for(Com in c("Root","Rhizosphere")){
IntDataSelectExp<- subset_samples(IntData, Experiment %in% Exp)
IntDataSelectExpCom<- subset_samples(IntDataSelectExp, Compartment %in% Com)
OtuTable<-data.frame(otu_table(IntDataSelectExpCom))
SampleData<-data.frame(sample_data(IntDataSelectExpCom))
library(vegan)
library(dplyr)
data = na.omit(OtuTable)%>%
t() %>%
as.data.frame() %>%
mutate(Variety = SampleData$Variety) %>%
mutate(Timepoint = SampleData$Timepoint) %>%
mutate(Treatment = SampleData$Treatment)  
site = data.frame(sample = rownames(t(OtuTable)),Variety = factor(as.character(data$Variety)),Timepoint = factor(as.character(data$Timepoint)),Treatment = factor(as.character(data$Treatment)))
ad<-adonis(formula = t(OtuTable) ~ Variety*Timepoint*Treatment, data = site,      permutations = 999, method = "bray")
#ad[[1]]
write.csv(ad[[1]],paste("OutputDir/Table1",Exp,Com,"Permanova.csv",sep="_"))
}
}



