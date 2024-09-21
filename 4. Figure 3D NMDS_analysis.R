rm(list = ls())
library(ggplot2)
library(phyloseq)
library(plyr)
mytheme3<-theme_bw()+
          theme(legend.position="right",
                             panel.border=element_rect(color="black", size=0.5, linetype="solid"),
                             panel.grid.major=element_line(linetype="dashed"),
                             panel.grid.minor=element_blank(),
                             plot.title=element_text(size=15,
                                                     colour="#003087",
                                                     family="CA"),
                             legend.text=element_text(size=9,colour="#003087",
                                                      family="CA"),
                             legend.key=element_blank(),
                             axis.text=element_text(size=10,colour="#003087",
                                                    family="CA"),
                             strip.text=element_text(size=12,colour="#EF0808",
                                                     family="CA"),
                             strip.background=element_blank()
                             
                )
mypoint2=geom_point(size=2)

IntData = readRDS("IntegratedBacteriaData.rds")
set.seed(123)##设计一个随机种子便于重复
IntDataRarefy<-rarefy_even_depth(IntData,rngseed=1,sample.size=15033,replace=T)
SampleData<-data.frame(sample_data(IntDataRarefy))
SampleData$Experiment<-factor(SampleData$Experiment,level=c("Exp1","Exp2"))
SampleData$Compartment<-factor(SampleData$Compartment,level=c("Root","Rhizosphere","BS"))
SampleData$MergeComExp<-paste(SampleData$Compartment,SampleData$Experiment,sep="_")
sample_data(IntDataRarefy)<-SampleData

GP.NMDS <- ordinate(IntDataRarefy, "NMDS")
Plot <- plot_ordination(IntDataRarefy,  GP.NMDS  , "samples", color="MergeComExp") +    geom_point(size=1.5)+ theme_classic()+  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

mydata <- data.frame(MergeComExp =Plot$data$MergeComExp,dn =Plot$data$NMDS1,dc =Plot$data$NMDS2)
chulls <- ddply(mydata, .(MergeComExp), function(mydata) mydata[chull(mydata$dn, mydata$dc), ])
scaleFUN <- function(x) sprintf("%.1f", x)###坐标轴小数点位数

FinalPlot<-Plot+mytheme3+mypoint2+ 
     scale_y_continuous(labels=scaleFUN)+ 
	 scale_x_continuous(labels=scaleFUN)+   
	 geom_polygon(data=chulls, aes(x=dn, y=dc, fill=MergeComExp, alpha=0.1 ))+ 
	 geom_point(aes( shape= Experiment),size=4)+
     scale_colour_manual(values = c("Root_Exp1"='SeaGreen4',"Root_Exp2"='SeaGreen4',"Rhizosphere_Exp1"='DarkSeaGreen3',"Rhizosphere_Exp2"='DarkSeaGreen3',"BS_Exp1"='Grey',"BS_Exp2"='Grey'))+
     scale_fill_manual(values = c("Root_Exp1"='SeaGreen4',"Root_Exp2"='SeaGreen4',"Rhizosphere_Exp1"='DarkSeaGreen3',"Rhizosphere_Exp2"='DarkSeaGreen3',"BS_Exp1"='Grey',"BS_Exp2"='Grey'))+
     theme_classic()+
     theme(plot.title =element_text(hjust=0.5,size=24,face = "bold"),axis.title.x =element_text(color="black",size=22.5,face = "bold"), axis.title.y=element_text(color="black",size=22.5,face = "bold"),axis.text=element_text(color="black",size=22.5,face = "bold"))+
     theme(axis.line.x=element_line(color="black",size=1.3,lineend = 10),#x轴线
          axis.ticks.x=element_line(color="black",size=1.3,lineend = 10),#x刻度
          axis.line.y=element_line(color="black",size=1.3,lineend = 10),#y轴线
          axis.ticks.y=element_line(color="black",size=1.3,lineend = 10))+
		  theme(legend.position = "none")


pdf("OutputDir/FinalPlot_NMDS.pdf")
FinalPlot
dev.off()

