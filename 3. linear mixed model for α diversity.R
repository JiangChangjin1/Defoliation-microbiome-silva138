
set.seed(5)
IntData = readRDS("IntegratedBacteriaData.rds")

erie<-IntData
min_lib <- min(sample_sums(erie))

nsamp=nsamples(erie)

trials=100

richness<-matrix(nrow=nsamp,ncol=trials)

row.names(richness)<-sample_names(erie)

evenness<-matrix(nrow=nsamp,ncol=trials)

row.names(evenness)<-sample_names(erie)

Shannon<-matrix(nrow=nsamp,ncol=trials)

row.names(Shannon)<-sample_names(erie)

ACE<-matrix(nrow=nsamp,ncol=trials)

row.names(ACE)<-sample_names(erie)

Fisher<-matrix(nrow=nsamp,ncol=trials)

row.names(Fisher)<-sample_names(erie)

Chao1<-matrix(nrow=nsamp,ncol=trials)

row.names(Chao1)<-sample_names(erie)

Simpson<-matrix(nrow=nsamp,ncol=trials)

row.names(Simpson)<-sample_names(erie)






set.seed(4)
for(i in 1:100){
#Subsample
r<-rarefy_even_depth(erie,
sample.size=min_lib,
verbose=FALSE,
replace=TRUE)

rich<-as.numeric(as.matrix(estimate_richness(r,measures='Observed')))

richness[,i]<-rich

even<-as.numeric(as.matrix(estimate_richness(r,measures='InvSimpson')))

evenness[,i]<-even

Shan<-as.numeric(as.matrix(estimate_richness(r,measures='Shannon')))

Shannon[,i]<-Shan

Cha<-as.numeric(as.matrix(estimate_richness(r,measures='Chao1')[,1]))

Chao1[,i]<-Cha

Sim<-as.numeric(as.matrix(estimate_richness(r,measures='Simpson')))

Simpson[,i]<-Sim


Fis<-as.numeric(as.matrix(estimate_richness(r,measures='Fisher')))

Fisher[,i]<-Fis

A<-as.numeric(as.matrix(estimate_richness(r,measures='ACE')[,1]))

ACE[,i]<-A

}
SampleID<-row.names(richness)
mean<-apply(richness,1,mean)
sd<-apply(richness,1,sd)
measure<-rep('Richness',nsamp)
rich_stats<-data.frame(SampleID,mean,sd,measure)


SampleID<-row.names(evenness)
mean<-apply(evenness,1,mean)
sd<-apply(evenness,1,sd)
measure<-rep('Inverse Simpson',nsamp)
even_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Shannon)
mean<-apply(Shannon,1,mean)
sd<-apply(Shannon,1,sd)
measure<-rep('Shannon',nsamp)
Shannon_stats<-data.frame(SampleID,mean,sd,measure)


SampleID<-row.names(ACE)
mean<-apply(ACE,1,mean)
sd<-apply(ACE,1,sd)
measure<-rep('ACE',nsamp)
ACE_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Simpson)
mean<-apply(Simpson,1,mean)
sd<-apply(Simpson,1,sd)
measure<-rep('Simpson',nsamp)
Simpson_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Fisher)
mean<-apply(Fisher,1,mean)
sd<-apply(Fisher,1,sd)
measure<-rep('Fisher',nsamp)
Fisher_stats<-data.frame(SampleID,mean,sd,measure)

SampleID<-row.names(Chao1)
mean<-apply(Chao1,1,mean)
sd<-apply(Chao1,1,sd)
measure<-rep('Chao1',nsamp)
Chao1_stats<-data.frame(SampleID,mean,sd,measure)

alphadiv<-cbind(as.character(rich_stats$SampleID),rich_stats$mean,even_stats$mean,Shannon_stats$mean,ACE_stats$mean,Simpson_stats$mean,Fisher_stats$mean,Chao1_stats$mean)

colnames(alphadiv)<-c("SampleID","observed_otu","evenness","Shannon","ACE","Simpson","Fisher","Chao1")
write.csv(alphadiv,"OutputDir/alphadiv.csv")








rm(list = ls())
 
library(nlme)                    # Fit Gaussian linear and nonlinear mixed-effects models
library(lme4)                    # Fit linear and generalized linear mixed-effects models
library(lattice)                 # Data visualization system
 

alphadiv<-read.csv("OutputDir/alphadiv.csv")

library(ggplot2)
library(phyloseq)
library(ggsignif)
scaleFUN <- function(x) sprintf("%.1f", x)
IntData = readRDS("IntegratedBacteriaData.rds")
AllLmeRelTable<-data.frame()
for( Exp in c("Exp1","Exp2")){
	IntDataExp<- subset_samples(IntData, Experiment %in% Exp)
	IntDataExpCom<- subset_samples(IntDataExp, Compartment %in% c("Root","Rhizosphere"))
	meta<-sample_data(IntDataExpCom)
	meta<-data.frame(meta[which(meta$SampleID %in% alphadiv$SampleID),])
	table<-merge(meta,alphadiv,by.x="SampleID")
		for(Val in c("Chao1","Simpson","Shannon","observed_otu")){
			if(Val == "Chao1"){
			LmeRel<-lme(Chao1~factor(Treatment)+factor(Variety)+factor(Compartment),random=list(Timepoint=~1),method="ML",data=table)
			LmeRelTable<-data.frame(summary(LmeRel)[20])
			LmeRelTable$group<-paste(Exp,Val,sep="_")
			AllLmeRelTable<-rbind(AllLmeRelTable,LmeRelTable)
			errorbar_up<-function(x){
                      mean(x)+sd(x)
            }
            errorbar_down<-function(x){
                      mean(x)-sd(x)
            }
			table$Treatment<-factor(table$Treatment,level=c("Intact","Defoliation"))
			p0<-ggplot(table,aes(x=Treatment,y=Chao1,color=Treatment,fill=Treatment))+
				geom_bar(stat = 'summary',fun=mean,width=0.6)+
				stat_summary(geom="errorbar", 
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color="black")+
				facet_wrap(~Timepoint*Compartment*Variety,nrow=1)+
				theme_bw() + 
				ylab(paste(Val,"index",sep=" "))+
				theme(panel.grid=element_blank(),
				legend.position = "none",
				axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
				coord_cartesian(clip = "off")+
				scale_color_manual(values =c("#406CB4","#E37933"))+
				scale_fill_manual(values =c("#406CB4","#E37933"))+
				theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
				theme(axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))#+
				ggsave(paste("OutputDir/lmeFigureS5/LME",Val,Exp,"alpha_index_bar.pdf",sep="_"),p0,width=15,height=5)
			}

			if(Val == "Simpson"){
			LmeRel<-lme(Simpson~factor(Treatment)+factor(Variety)+factor(Compartment),random=list(Timepoint=~1),method="ML",data=table)
			LmeRelTable<-data.frame(summary(LmeRel)[20])
			LmeRelTable$group<-paste(Exp,Val,sep="_")
			AllLmeRelTable<-rbind(AllLmeRelTable,LmeRelTable)
			errorbar_up<-function(x){
                      mean(x)+sd(x)
            }
            errorbar_down<-function(x){
                      mean(x)-sd(x)
            }
			table$Treatment<-factor(table$Treatment,level=c("Intact","Defoliation"))
			p0<-ggplot(table,aes(x=Treatment,y=Simpson,color=Treatment,fill=Treatment))+
				geom_bar(stat = 'summary',fun=mean,width=0.6)+
				stat_summary(geom="errorbar", 
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color="black")+
				facet_wrap(~Timepoint*Compartment*Variety,nrow=1)+
				theme_bw() + 
				ylab(paste(Val,"index",sep=" "))+
				theme(panel.grid=element_blank(),
				legend.position = "none",
				axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
				coord_cartesian(clip = "off")+
				scale_color_manual(values =c("#406CB4","#E37933"))+
				scale_fill_manual(values =c("#406CB4","#E37933"))+
				theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
				theme(axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))#+
				ggsave(paste("OutputDir/lmeFigureS5/LME",Val,Exp,"alpha_index_bar.pdf",sep="_"),p0,width=15,height=5)
			}
			if(Val == "Shannon"){
			LmeRel<-lme(Shannon~factor(Treatment)+factor(Variety)+factor(Compartment),random=list(Timepoint=~1),method="ML",data=table)
			LmeRelTable<-data.frame(summary(LmeRel)[20])
			LmeRelTable$group<-paste(Exp,Val,sep="_")
			AllLmeRelTable<-rbind(AllLmeRelTable,LmeRelTable)
			errorbar_up<-function(x){
                      mean(x)+sd(x)
            }
            errorbar_down<-function(x){
                      mean(x)-sd(x)
            }
			table$Treatment<-factor(table$Treatment,level=c("Intact","Defoliation"))
			p0<-ggplot(table,aes(x=Treatment,y=Shannon,color=Treatment,fill=Treatment))+
				geom_bar(stat = 'summary',fun=mean,width=0.6)+
				stat_summary(geom="errorbar", 
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color="black")+
				facet_wrap(~Timepoint*Compartment*Variety,nrow=1)+
				theme_bw() + 
				ylab(paste(Val,"index",sep=" "))+
				theme(panel.grid=element_blank(),
				legend.position = "none",
				axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
				coord_cartesian(clip = "off")+
				scale_color_manual(values =c("#406CB4","#E37933"))+
				scale_fill_manual(values =c("#406CB4","#E37933"))+
				theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
				theme(axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))#+
				ggsave(paste("OutputDir/lmeFigureS5/LME",Val,Exp,"alpha_index_bar.pdf",sep="_"),p0,width=15,height=5)
			}
			if(Val == "observed_otu"){
			LmeRel<-lme(observed_otu~factor(Treatment)+factor(Variety)+factor(Compartment),random=list(Timepoint=~1),method="ML",data=table)
			LmeRelTable<-data.frame(summary(LmeRel)[20])
			LmeRelTable$group<-paste(Exp,Val,sep="_")
			AllLmeRelTable<-rbind(AllLmeRelTable,LmeRelTable)
			errorbar_up<-function(x){
                      mean(x)+sd(x)
            }
            errorbar_down<-function(x){
                      mean(x)-sd(x)
            }
			table$Treatment<-factor(table$Treatment,level=c("Intact","Defoliation"))
			p0<-ggplot(table,aes(x=Treatment,y=observed_otu,color=Treatment,fill=Treatment))+
				geom_bar(stat = 'summary',fun=mean,width=0.6)+
				stat_summary(geom="errorbar", 
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color="black")+
				facet_wrap(~Timepoint*Compartment*Variety,nrow=1)+
				theme_bw() + 
				ylab(paste(Val,"index",sep=" "))+
				theme(panel.grid=element_blank(),
				legend.position = "none",
				axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
				coord_cartesian(clip = "off")+
				scale_color_manual(values =c("#406CB4","#E37933"))+
				scale_fill_manual(values =c("#406CB4","#E37933"))+
				theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
				theme(axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))#+
				ggsave(paste("OutputDir/lmeFigureS5/LME",Val,Exp,"alpha_index_bar.pdf",sep="_"),p0,width=15,height=5)
			}
		}
	}


write.csv(AllLmeRelTable,"OutputDir/lmeFigureS5/LME_alpha_index.csv")









