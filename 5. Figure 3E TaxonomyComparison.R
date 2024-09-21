rm(list = ls())
library("rlang")
library(ggplot2)
library(phyloseq)
scaleFUN <- function(x) sprintf("%.1f", x)
IntData = readRDS("IntegratedBacteriaData.rds")
IntData = transform_sample_counts(IntData, function(x) x / sum(x) )
IntData<- subset_samples(IntData, !Treatment %in% c("Defoliation"))

OtuTable<-data.frame(otu_table(IntData))
TaxTable<-data.frame(tax_table(IntData))
OtuTable$id<-rownames(OtuTable)
TaxTable$id<- rownames(TaxTable)
MergeOtuTaxTable = merge(OtuTable,TaxTable, by = "id")
MergeOtuTaxTable$Phylum<-as.vector(MergeOtuTaxTable$Phylum)
MergeOtuTaxTable$Class<-as.vector(MergeOtuTaxTable$Class)
MergeOtuTaxTable[which(MergeOtuTaxTable$Phylum == "Proteobacteria"),]$Phylum<-MergeOtuTaxTable[which(MergeOtuTaxTable$Phylum == "Proteobacteria"),]$Class
num<-which(colnames(MergeOtuTaxTable)=="Phylum")
PhySumTable<-aggregate(MergeOtuTaxTable[,2:(which(colnames(MergeOtuTaxTable)=="Kingdom")-1)], by=MergeOtuTaxTable[num], FUN=sum)
rownames(PhySumTable)<-PhySumTable[,1]
PhySumTable<-PhySumTable[,-1]
PhySumTable$sum <- rowSums(PhySumTable)
#求各类群的丰度总和，并排序
PhySumTable <-  PhySumTable[order(PhySumTable$sum,decreasing = TRUE),]
#删除最后一列和,并取Top10物种作图
PhySumTableTop10 <- PhySumTable[c(1:10),-ncol(PhySumTable)]
# 剩余的物种合并为Others
PhySumTableTop10['Others',] <- 1-colSums(PhySumTableTop10)
PhySumTableTop10<-data.frame(t(PhySumTableTop10))
SampleData<-data.frame(sample_data(IntData))
SampleData$MergeExpCom<-paste(SampleData$Experiment,SampleData$Compartment,sep="_")
PhySumTableTop10$SampleID<-rownames(PhySumTableTop10)
PhySumTableTop10MergeInfo = merge(PhySumTableTop10,SampleData, by = "SampleID")
for(i in 2:(which(colnames(PhySumTableTop10MergeInfo)=="Experiment")-1)){
  PhySumTableTop10MergeInfo[,i]<-as.numeric(as.vector(PhySumTableTop10MergeInfo[,i]))
}



num<-which(colnames(PhySumTableTop10MergeInfo)=="MergeExpCom")
PhySumTableTop10MergeInfoMean<-aggregate(PhySumTableTop10MergeInfo[,2:(which(colnames(PhySumTableTop10MergeInfo)=="Experiment")-1)], by=PhySumTableTop10MergeInfo[num], FUN=mean)
rownames(PhySumTableTop10MergeInfoMean)<-PhySumTableTop10MergeInfoMean[,1]
PhySumTableTop10MergeInfoMean<-PhySumTableTop10MergeInfoMean[,-1]
PhySumTableTop10MergeInfoMean$sample<-rownames(PhySumTableTop10MergeInfoMean)
PhySumTableTop10MergeInfoMean<-PhySumTableTop10MergeInfoMean[c(1,4,2,5,3,6),]
library(reshape2)
 SampleTaxAbundance<- melt(PhySumTableTop10MergeInfoMean, id.vars = c('sample'), variable.name = 'Phylum', value.name = 'abundance')
unique(SampleTaxAbundance$Phylum)
 [1] Gammaproteobacteria Myxococcota         Bacteroidota        Desulfobacterota    Alphaproteobacteria Firmicutes          Spirochaetota      
 [8] Nitrospirota        Verrucomicrobiota   Acidobacteriota     Others             
11 Levels: Gammaproteobacteria Myxococcota Bacteroidota Desulfobacterota Alphaproteobacteria Firmicutes Spirochaetota ... Others

#result2 <- melt(phylum_top10_t3_2, id.vars = c('sample'), variable.name = 'Phylum', value.name = 'abundance')
SampleTaxAbundance$Phylum<-factor(SampleTaxAbundance$Phylum,levels=c("Gammaproteobacteria","Myxococcota","Bacteroidota","Desulfobacterota","Alphaproteobacteria","Firmicutes","Spirochaetota","Nitrospirota","Verrucomicrobiota","Acidobacteriota","Others"))
SampleTaxAbundance$sample<-factor(SampleTaxAbundance$sample,levels=unique(PhySumTableTop10MergeInfoMean$sample))


Plot <-ggplot(SampleTaxAbundance,
              aes(x=sample,
                  y=abundance*100,
                  fill=Phylum,
                  stratum = Phylum,
                  alluvium = Phylum)) +
  geom_bar(stat='identity', width=0.45) +
  #geom_alluvium() +
  #geom_stratum(width=0.45, size=0.5) +
  labs(x='Samples', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_y_continuous(labels=scaleFUN)+
  theme(plot.title =element_text(hjust=0.5,size=24,color="black",face = "bold"),axis.title.x =element_text(size=22.5,color="black",face = "bold"), axis.title.y=element_text(size=22.5,color="black",face = "bold"),axis.text=element_text(size=22.5,color="black",face = "bold"))+
  theme(axis.line.x=element_blank(),#x轴线
        axis.ticks.x=element_blank(),#x刻度
        axis.line.y=element_line(color="black",size=1.3,lineend = 10),#y轴线
        axis.ticks.y=element_line(color="black",size=1.3,lineend = 10),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank())+scale_y_continuous(breaks = c(0,100),expand = c(0,0))+
  scale_fill_manual(values = c(Gammaproteobacteria='RoyalBlue',Myxococcota='Burlywood3',Bacteroidota='IndianRed1',
                               Desulfobacterota='#9999FF',Alphaproteobacteria='DarkGoldenrod4',Firmicutes='MistyRose3',Spirochaetota='#d1bbff',Nitrospirota='#770077',Verrucomicrobiota='#ee7700',
                               Acidobacteriota='#cceeff',Others="gray")
  )

pdf("OutputDir/TaxonomyComparison.pdf")
Plot
dev.off()

