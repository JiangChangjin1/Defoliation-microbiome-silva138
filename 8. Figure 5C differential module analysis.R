

###Rhizosphere
rm(list = ls())
set.seed(5)
library(phyloseq)
library(ggpubr)

IntData = readRDS("IntegratedBacteriaData.rds")
tab_all<-c()
tab_all2<-c()
rel_all_module<-c()
for(Com in c("Rhizosphere")){
	for(Exp in c("Exp1","Exp2")){
		IntDataSelCom<- subset_samples(IntData, Compartment %in% Com)
		IntDataSelComExp<- subset_samples(IntDataSelCom, Experiment %in% Exp)
		IntDataSelComExp <- transform_sample_counts(IntDataSelComExp, function(x) x / sum(x) )
		for(Var in c("NIP","MH63")){
			IntDataSelComExpVar<- subset_samples(IntDataSelComExp, Variety %in% Var )
			IntDataSelComExpVar<- subset_samples(IntDataSelComExpVar, Treatment %in% c("Intact","Defoliation"))
			sample_data<-data.frame(sample_data(IntDataSelComExpVar))
			sample_data$Treatment <- factor(sample_data$Treatment,levels=c("Intact","Defoliation"))
			sample_data$Timepoint <- factor(sample_data$Timepoint,unique(sample_data$Timepoint))
			otu_table<-data.frame(otu_table(IntDataSelComExpVar))
			NodesTable<-read.csv("OutputDir/NetworkAnalysis/Add_module_Info/Nodes_Table_Rhizosphere_network_add_module_Info.csv")
			ModuleAbun<-c()#otu_table_m<-c()
			for(m in unique(NodesTable$modularity_class)){
			#rs_node1<-subset(NodesTable,modularity_class == m)
			NodesTableSel<-subset(NodesTable,modularity_class == m)
			otu_tableSelMod<-otu_table[which(rownames(otu_table) %in% NodesTableSel$Id),]
			if(length(which(colSums(otu_tableSelMod)>0))/ncol(otu_tableSelMod)>0.1){ ##Remove low-occurrence modules
				if(length(which(rowSums(otu_tableSelMod)==0))>0){
				otu_tableSelMod<-otu_tableSelMod[-which(rowSums(otu_tableSelMod)==0),]
				}
				zmodu1<-t(scale(t(otu_tableSelMod)))
				SelModuleAbun<-as.data.frame(colSums(zmodu1)/nrow(zmodu1))
				colnames(SelModuleAbun)<-m
				SelModuleAbunFinal<-t(SelModuleAbun)
				ModuleAbun<-rbind(ModuleAbun,SelModuleAbunFinal)
				}
			}
			colnames(ModuleAbun)<-colnames(otu_table)
			rownames(ModuleAbun)<-paste("module",rownames(ModuleAbun),sep="")
			count =  data.frame(t(ModuleAbun))
			count$SampleID<-rownames(count)
			table<-merge(count,sample_data,by.x="SampleID")
			rel_all<-c()
			for(n in 2:(which(colnames(table) == "Experiment")-1)){
			table0<-table[,c(n,which(colnames(table) == "Experiment"):ncol(table))]
			colnames(table0)[1]<-"module"
			p_all<-c()
			t_all<-c()
			def_mean<-c()
			int_mean<-c()
			for(t0 in unique(table0$Timepoint) ){
				unit.lm <- t.test(module  ~ Treatment , data = subset(table0, Timepoint  == t0 ) )
				p_all<-c(p_all,unit.lm$p.value)
				t_all<-c(t_all, unit.lm$statistic)
				def_mean<-c(def_mean,mean(subset(subset(table0, Timepoint  == t0 ),Treatment=="Defoliation")$module ))
				int_mean<-c(int_mean,mean(subset(subset(table0, Timepoint  == t0 ),Treatment=="Intact")$module ))
				}
				rel<-data.frame(timepoint=unique(table0$Timepoint),Pvalue=p_all,Tvalue=t_all,def_mean=def_mean,int_mean=int_mean)
				rel$module<-colnames(table)[n]
				rel_all<-rbind(rel_all,rel)
				}
			rel_all<-data.frame(rel_all)
			rel_all$Experiment<-Exp
			rel_all$Variety<-Var
			rel_all$Compartment<-Com
			rel_all_module<-rbind(rel_all_module,rel_all)
			}
		}
	}


write.csv(rel_all_module,"OutputDir/NetworkAnalysis/Add_module_Info/Rhizosphere_module_zscore_t.test_result.csv")
	  



rel_all_fun_sig<-subset(rel_all_module,Pvalue<0.05)
 unique(rel_all_fun_sig$module)
 [1] "module1"  "module13" "module8"  "module2"  "module3"  "module14" "module11" "module15" "module16"
[10] "module17" "module32" "module5"  "module6"  "module9"  "module20" "module22" "module24" "module26"
[19] "module30" "module10" "module21" "module7"  "module12" "module18" "module19" "module27" "module28"
[28] "module33"

all_tab<-c()
for(Com in c("Rhizosphere")){
    for(Var in c("NIP","MH63")){
		rel_all_fun_sig1<-subset(rel_all_fun_sig,Variety == Var)
		rel_all_fun_sig2<-subset(rel_all_fun_sig1,Compartment == Com)
		sel<-table(rel_all_fun_sig2$module )[which(table(rel_all_fun_sig2$module) > 1)]
		rel_all_fun_sig3<-rel_all_fun_sig2[which(rel_all_fun_sig2$module  %in% names(sel)),]
		rel_all_fun_sig_final<-c()
		for(f in unique(rel_all_fun_sig3$module )){
			rel_all_fun_sig4<-subset(rel_all_fun_sig3,module  == f)
			if(length(unique(rel_all_fun_sig4$Experiment))>1){
				rel_all_fun_sig_final<-rbind(rel_all_fun_sig_final,rel_all_fun_sig4)
				}
			}
		all_tab<-rbind(all_tab,rel_all_fun_sig_final)
		}
	}


all_tab_sel<-c()
for(Com in c("Rhizosphere")){
    for(Var in c("NIP","MH63")){
		all_tab01<-subset(all_tab,Compartment == Com)
		all_tab02<-subset(all_tab01,Variety == Var)
		for(f in unique(all_tab02$module )){
			all_tab1<-subset(all_tab02,module  == f)
			if((length(unique(all_tab1[which(all_tab1$Tvalue > 0),]$Experiment))>1)|(length(unique(all_tab1[which(all_tab1$Tvalue < 0),]$Experiment))>1)){
				all_tab_sel<-rbind(all_tab_sel,all_tab1)
				}
			}
		}
	}
	
all_tab_sel0<-c()
for(Var in c("MH63","NIP")){
	all_tab_sel_var<-subset(all_tab_sel,Variety == Var)
	for(i in unique(all_tab_sel$module)){
		if((length(which(subset(all_tab_sel_var,module == i)$Tvalue>0))==0)|(length(which(subset(all_tab_sel_var,module == i)$Tvalue<0))==0)){
			all_tab_sel0<-rbind(all_tab_sel0,subset(all_tab_sel_var,module == i))
			}
		}	
	}


write.csv(all_tab_sel0,"OutputDir/NetworkAnalysis/Add_module_Info/Rhizosphere_module_zscore_t.test_result_SelConExp.csv")











###Root
rm(list = ls())
set.seed(5)
library(phyloseq)
library(ggpubr)

IntData = readRDS("IntegratedBacteriaData.rds")
tab_all<-c()
tab_all2<-c()
rel_all_module<-c()
for(Com in c("Root")){
	for(Exp in c("Exp1","Exp2")){
		IntDataSelCom<- subset_samples(IntData, Compartment %in% Com)
		IntDataSelComExp<- subset_samples(IntDataSelCom, Experiment %in% Exp)
		IntDataSelComExp <- transform_sample_counts(IntDataSelComExp, function(x) x / sum(x) )
		for(Var in c("NIP","MH63")){
			IntDataSelComExpVar<- subset_samples(IntDataSelComExp, Variety %in% Var )
			IntDataSelComExpVar<- subset_samples(IntDataSelComExpVar, Treatment %in% c("Intact","Defoliation"))
			sample_data<-data.frame(sample_data(IntDataSelComExpVar))
			sample_data$Treatment <- factor(sample_data$Treatment,levels=c("Intact","Defoliation"))
			sample_data$Timepoint <- factor(sample_data$Timepoint,unique(sample_data$Timepoint))
			otu_table<-data.frame(otu_table(IntDataSelComExpVar))
			NodesTable<-read.csv("OutputDir/NetworkAnalysis/Add_module_Info/Nodes_Table_Root_network_add_module_Info.csv")
			ModuleAbun<-c()#otu_table_m<-c()
			for(m in unique(NodesTable$modularity_class)){
			#rs_node1<-subset(NodesTable,modularity_class == m)
			NodesTableSel<-subset(NodesTable,modularity_class == m)
			otu_tableSelMod<-otu_table[which(rownames(otu_table) %in% NodesTableSel$Id),]
			if(length(which(colSums(otu_tableSelMod)>0))/ncol(otu_tableSelMod)>0.1){ ##Remove low-occurrence modules
				if(length(which(rowSums(otu_tableSelMod)==0))>0){
				otu_tableSelMod<-otu_tableSelMod[-which(rowSums(otu_tableSelMod)==0),]
				}
				zmodu1<-t(scale(t(otu_tableSelMod)))
				SelModuleAbun<-as.data.frame(colSums(zmodu1)/nrow(zmodu1))
				colnames(SelModuleAbun)<-m
				SelModuleAbunFinal<-t(SelModuleAbun)
				ModuleAbun<-rbind(ModuleAbun,SelModuleAbunFinal)
				}
			}
			colnames(ModuleAbun)<-colnames(otu_table)
			rownames(ModuleAbun)<-paste("module",rownames(ModuleAbun),sep="")
			count =  data.frame(t(ModuleAbun))
			count$SampleID<-rownames(count)
			table<-merge(count,sample_data,by.x="SampleID")
			rel_all<-c()
			for(n in 2:(which(colnames(table) == "Experiment")-1)){
			table0<-table[,c(n,which(colnames(table) == "Experiment"):ncol(table))]
			colnames(table0)[1]<-"module"
			p_all<-c()
			t_all<-c()
			def_mean<-c()
			int_mean<-c()
			for(t0 in unique(table0$Timepoint) ){
				unit.lm <- t.test(module  ~ Treatment , data = subset(table0, Timepoint  == t0 ) )
				p_all<-c(p_all,unit.lm$p.value)
				t_all<-c(t_all, unit.lm$statistic)
				def_mean<-c(def_mean,mean(subset(subset(table0, Timepoint  == t0 ),Treatment=="Defoliation")$module ))
				int_mean<-c(int_mean,mean(subset(subset(table0, Timepoint  == t0 ),Treatment=="Intact")$module ))
				}
				rel<-data.frame(timepoint=unique(table0$Timepoint),Pvalue=p_all,Tvalue=t_all,def_mean=def_mean,int_mean=int_mean)
				rel$module<-colnames(table)[n]
				rel_all<-rbind(rel_all,rel)
				}
			rel_all<-data.frame(rel_all)
			rel_all$Experiment<-Exp
			rel_all$Variety<-Var
			rel_all$Compartment<-Com
			rel_all_module<-rbind(rel_all_module,rel_all)
			}
		}
	}


write.csv(rel_all_module,"OutputDir/NetworkAnalysis/Add_module_Info/Root_module_zscore_t.test_result.csv")
	  



rel_all_fun_sig<-subset(rel_all_module,Pvalue<0.05)
[1] "module8"  "module6"  "module17" "module13" "module15" "module1"  "module18" "module9"  "module11"
[10] "module14" "module20" "module16" "module10" "module19"


all_tab<-c()
for(Com in c("Root")){
    for(Var in c("NIP","MH63")){
		rel_all_fun_sig1<-subset(rel_all_fun_sig,Variety == Var)
		rel_all_fun_sig2<-subset(rel_all_fun_sig1,Compartment == Com)
		sel<-table(rel_all_fun_sig2$module )[which(table(rel_all_fun_sig2$module) > 1)]
		rel_all_fun_sig3<-rel_all_fun_sig2[which(rel_all_fun_sig2$module  %in% names(sel)),]
		rel_all_fun_sig_final<-c()
		for(f in unique(rel_all_fun_sig3$module )){
			rel_all_fun_sig4<-subset(rel_all_fun_sig3,module  == f)
			if(length(unique(rel_all_fun_sig4$Experiment))>1){
				rel_all_fun_sig_final<-rbind(rel_all_fun_sig_final,rel_all_fun_sig4)
				}
			}
		all_tab<-rbind(all_tab,rel_all_fun_sig_final)
		}
	}


all_tab_sel<-c()
for(Com in c("Root")){
    for(Var in c("NIP","MH63")){
		all_tab01<-subset(all_tab,Compartment == Com)
		all_tab02<-subset(all_tab01,Variety == Var)
		for(f in unique(all_tab02$module )){
			all_tab1<-subset(all_tab02,module  == f)
			if((length(unique(all_tab1[which(all_tab1$Tvalue > 0),]$Experiment))>1)|(length(unique(all_tab1[which(all_tab1$Tvalue < 0),]$Experiment))>1)){
				all_tab_sel<-rbind(all_tab_sel,all_tab1)
				}
			}
		}
	}
	
all_tab_sel0<-c()
for(Var in c("MH63","NIP")){
	all_tab_sel_var<-subset(all_tab_sel,Variety == Var)
	for(i in unique(all_tab_sel$module)){
		if((length(which(subset(all_tab_sel_var,module == i)$Tvalue>0))==0)|(length(which(subset(all_tab_sel_var,module == i)$Tvalue<0))==0)){
			all_tab_sel0<-rbind(all_tab_sel0,subset(all_tab_sel_var,module == i))
			}
		}	
	}


write.csv(all_tab_sel0,"OutputDir/NetworkAnalysis/Add_module_Info/Root_module_zscore_t.test_result_SelConExp.csv")















