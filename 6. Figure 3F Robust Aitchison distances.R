		
			
			rm(list = ls())
            beta_data_all_tre<-c()
            library(phyloseq)
            library(ggsignif)
			library(ggplot2)
			library(tidyr)
            
            
            IntData = readRDS("IntegratedBacteriaData.rds")
            IntDataRar<-rarefy_even_depth(IntData,rngseed=1,sample.size= min(colSums(otu_table(IntData))),replace=F)
            IntDataRar<- subset_samples(IntDataRar, !Compartment %in% "BS")

            beta_data_all_tre<-c()
            for(Exp in c("Exp1","Exp2")){
              for(Com in c("Root","Rhizosphere")){
                for(Var in c("NIP","MH63")){
                  
                  if(Exp == "Exp1"){
                    time_all<-c("day1","day3","day5","day7")
                  }
                  if(Exp == "Exp2"){
                    time_all<-c("day2","day3","day4","day5")  
                  }
                  
                  for(Tim in time_all){
                    
                    IntDataRarExp<- subset_samples(IntDataRar, Experiment %in% Exp)
                    IntDataRarExpCom<- subset_samples(IntDataRarExp, Compartment %in% Com)
                    IntDataRarExpComVar<- subset_samples(IntDataRarExpCom, Variety %in% Var)
                    IntDataRarExpComVarTim<- subset_samples(IntDataRarExpComVar, Timepoint %in% Tim)
                    otu<-otu_table(IntDataRarExpComVarTim)
                    otu_table<-t(otu)
                    bray_curtis = vegan::vegdist(otu_table, method = "robust.aitchison")
                    bray_curtis= as.matrix(bray_curtis)
                    IntDataRarExpComVarTim_1<- subset_samples(IntDataRarExpComVarTim, Treatment %in% "Intact")
                    s1<-rownames(sample_data(IntDataRarExpComVarTim_1))
                    IntDataRarExpComVarTim_2<- subset_samples(IntDataRarExpComVarTim, Treatment %in% "Defoliation")
                    s2<-rownames(sample_data(IntDataRarExpComVarTim_2))
                    g1<-unique(as.vector(bray_curtis[s1,s1]))
                    g1<-g1[-which(g1 == 0)]
                    g3<-unique(as.vector(bray_curtis[s1,s2]))
                    beta_data0<-data.frame(distance=c(g1,g3),group=paste(Exp,Com,Var,Tim,sep="_"),beta_group=rep(c("Within treatments","Between treatments"),c(length(g1) ,length(g3) )))
                    
                    beta_data_all_tre<-rbind(beta_data_all_tre,beta_data0)
                    
                  }
                }
              }
            }
            
            

            

            
		
            
            beta_data_all_tre<- separate(beta_data_all_tre, group, into = c("Experiment", "Compartment", "Variety","Timepoint"), sep = "_") 

            beta_data_all_tre$group_comb<-paste(beta_data_all_tre$Timepoint,beta_data_all_tre$beta_group,sep=" ")
            
            
            
            p<-list()
            i=1
            for(Exp in c("Exp1","Exp2")){ 
              for(Var in c("NIP","MH63")){ 
                for(Com in c("Root","Rhizosphere")){ 
                  beta_data_all_treExp<- subset(beta_data_all_tre, Experiment == Exp)
                  beta_data_all_treExpVar<- subset(beta_data_all_treExp, Variety == Var)
                  beta_data_all_treExpVarCom<- subset(beta_data_all_treExpVar, Compartment == Com)
                  if(Exp == "Exp1"){
                    beta_data_all_treExpVarCom$group_comb<-factor(beta_data_all_treExpVarCom$group_comb,levels = c("day1 Within treatments","day1 Between treatments","day3 Within treatments","day3 Between treatments","day5 Within treatments","day5 Between treatments","day7 Within treatments","day7 Between treatments")
                    )
                  }
                  
                  if(Exp == "Exp2"){
                    beta_data_all_treExpVarCom$group_comb<-factor(beta_data_all_treExpVarCom$group_comb,levels =  c("day2 Within treatments","day2 Between treatments","day3 Within treatments","day3 Between treatments","day4 Within treatments","day4 Between treatments","day5 Within treatments","day5 Between treatments")
                    )
                  }
                  
                  
                  
                  
                  if(Exp == "Exp1"){
                    errorbar_up<-function(x){ 
                      mean(x)+sd(x)
                    }
                    errorbar_down<-function(x){
                      mean(x)-sd(x)
                    }
                    p[[i]]<-ggplot(beta_data_all_treExpVarCom,aes(x = group_comb, y = distance, color = beta_group,fill=beta_group)) +
							geom_bar(stat = 'summary',width=0.6,fun=mean)+
							stat_summary(geom="errorbar",  
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color = "black")+
                      geom_signif(
                        textsize  = 4,
                        comparisons = list(
                          c("day1 Within treatments","day1 Between treatments"),
                          c("day3 Within treatments","day3 Between treatments"),
                          c("day5 Within treatments","day5 Between treatments"),
                          c("day7 Within treatments","day7 Between treatments")),
                        map_signif_level = T,
                        test = t.test, 
                        vjust=-0.2, 
						tip_length = 0,col="black", y_position =60 
                      )+
                      theme_bw() + 
                      xlab(paste(Exp,Var,Com,sep="_")) +
                      ylab("Robust Aitchison distance")+
                      theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+#+#+
                     ylim(0,65)+
                      coord_cartesian(clip = "off")+ 
                      scale_color_manual(values =c("#E37933","#406CB4"))+
                      scale_fill_manual(values =c("#E37933","#406CB4"))+
                          theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                      
                      theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))+
                      geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")#+
                    
                  }
                  
                  
                  
                  if(Exp == "Exp2"){
                    errorbar_up<-function(x){ 
                      mean(x)+sd(x)
                    }
                    errorbar_down<-function(x){
                      mean(x)-sd(x)
                    }
                    p[[i]]<-ggplot(beta_data_all_treExpVarCom,aes(x = group_comb, y = distance, color = beta_group,fill=beta_group)) +
                         geom_bar(stat = 'summary',width=0.6,fun=mean)+
                      stat_summary(geom="errorbar",   
                                   fun.min = errorbar_down,
                                   fun.max = errorbar_up,
                                   width=0.3,color = "black")+
                      geom_signif(
                        textsize  = 4,
                        comparisons = list(c("day2 Within treatments","day2 Between treatments"),
                                           c("day3 Within treatments","day3 Between treatments"),
                                           c("day4 Within treatments","day4 Between treatments"),
                                           c("day5 Within treatments","day5 Between treatments")),
                        map_signif_level = T,
                        test = t.test,
                        vjust=-0.2,
						tip_length = 0,col="black", y_position =60
                      )+
                      theme_bw() + 
                      xlab(paste(Exp,Var,Com,sep="_")) +
                      ylab("Robust Aitchison distance")+
                      theme(panel.grid=element_blank(),
                            legend.position = "none",
                            axis.text.x = element_text(angle=90, hjust=1, vjust=.6,size=6))+
                      ylim(0,65)+
                      coord_cartesian(clip = "off")+
                      scale_color_manual(values =c("#E37933","#406CB4"))+
                      scale_fill_manual(values =c("#E37933","#406CB4"))+
                      theme(plot.margin = unit(c(2,0.2,0.2,0.2),'cm'))+
                      
                      theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"))+
                      geom_vline(xintercept = seq(2.5,8,by=2),linetype="dotted")

                  }
                  
                  
                  
                  i=i+1
                  
                }
              }
            }
      
            
            library(Rmisc)
            pdf("OutputDir/RobustAitchisonDistances_Figure3E/robust.aitchison_modify_color.t.test.15x7.within.intact.pdf",width=10,height=5)
            multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],cols =4)
            
            dev.off()