
### ancombc2 pipline for Timepoint + Treatment

rm(list = ls())
set.seed(5)
library(phyloseq)
library(edgeR)
library(ANCOMBC)
IntData = readRDS("IntegratedBacteriaData.rds")
min_lib <- min(sample_sums(IntData))

for(Exp in c("Exp1","Exp2")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			IntDataSelExp<- subset_samples(IntData, Experiment %in% Exp)
			IntDataSelExpCom<- subset_samples(IntDataSelExp, Compartment %in% Com)
			IntDataSelExpCom<- subset_samples(IntDataSelExpCom, Treatment %in% c("Intact","Defoliation"))
			IntDataSelExpComVar<- subset_samples(IntDataSelExpCom, Variety %in% Var )
			tse = mia::makeTreeSummarizedExperimentFromPhyloseq(IntDataSelExpComVar)
			set.seed(123)
			out_Family = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
			fix_formula = "Timepoint + Treatment",
			rand_formula = NULL,
			p_adj_method = "fdr", pseudo_sens = TRUE,
			prv_cut = 0.10, lib_cut = min_lib, s0_perc = 0.05,
			group = "Treatment", struc_zero = TRUE, neg_lb = TRUE,
			alpha = 0.05, n_cl = 1, #verbose = TRUE,
			lme_control = lme4::lmerControl(),
			mdfdr_control = list(fwer_ctrl_method = "fdr", B = 1))#,
			saveRDS(object = out_Family,file = paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin.rds",sep="_"))
			}
		}
	}



### ancombc2 pipline for Timepoint * Treatment

rm(list = ls())
set.seed(5)
library(phyloseq)
library(edgeR)
library(ANCOMBC)
IntData  = readRDS("IntegratedBacteriaData.rds")
min_lib <- min(sample_sums(IntData))


for(Exp in c("Exp1","Exp2")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			IntDataSelExp<- subset_samples(IntData, Experiment %in% Exp)
			IntDataSelExpCom<- subset_samples(IntDataSelExp, Compartment %in% Com)
			IntDataSelExpCom<- subset_samples(IntDataSelExpCom, Treatment %in% c("Intact","Defoliation"))
			IntDataSelExpComVar<- subset_samples(IntDataSelExpCom, Variety %in% Var )
			tse = mia::makeTreeSummarizedExperimentFromPhyloseq(IntDataSelExpComVar)
			set.seed(123)
			out_Family = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
			fix_formula = "Timepoint * Treatment",
			rand_formula = NULL,
			p_adj_method = "fdr", pseudo_sens = TRUE,
			prv_cut = 0.10, lib_cut = min_lib, s0_perc = 0.05,
			group = "Treatment", struc_zero = TRUE, neg_lb = TRUE,
			alpha = 0.05, n_cl = 1, #verbose = TRUE,
			lme_control = lme4::lmerControl(),
			mdfdr_control = list(fwer_ctrl_method = "fdr", B = 1))#,
			saveRDS(object = out_Family,file = paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin_inter.rds",sep=""))

			}
		}
	}


## significant result for defoliation group
library(phyloseq)
IntData = readRDS("IntegratedBacteriaData.rds")
tab_all<-data.frame()
for(Exp in c("Exp1","Exp2")){
	for(Com in c("Root","Rhizosphere")){
		for(Var in c("MH63","NIP")){
			e1<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin.rds",sep=""))
			e1_sig<-e1$res[which(e1$res$p_TreatmentIntact < 0.01),]#$taxon
			e1_sig_n<-e1_sig[which(abs(e1_sig$lfc_TreatmentIntact) > log2(1.5)),]#

			tax<-data.frame(tax_table(IntData)[e1_sig_n$taxon,])
			tax$taxon<-rownames(tax)

			e1_sig_tax<-merge(e1_sig_n,tax, by.x = "taxon")

			e1_sig_tax$Experiment<-Exp
			e1_sig_tax$Compartment<-Com
			e1_sig_tax$Variety<-Var
			colnames(e1_sig_tax)[c(3:5,8:10,13:15,18:20,23:25,28:30)]<-   c("lfc_Timepoint2","lfc_Timepoint3","lfc_Timepoint4","se_Timepoint2","se_Timepoint3","se_Timepoint4","W_Timepoint2","W_Timepoint3","W_Timepoint4" ,  "p_Timepoint2","p_Timepoint3" ,  "p_Timepoint4" ,  "q_Timepoint2" ,  "q_Timepoint3" ,  "q_Timepoint4"  ,  "diff_Timepoint2" ,  "diff_Timepoint3" ,  "diff_Timepoint4") 
			tab_all<-rbind(tab_all,e1_sig_tax)
			}
		}
	}

write.csv(tab_all,"OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_138.csv")




### significant result for sample timeponit group
rm(list = ls())
library(phyloseq)
IntData = readRDS("IntegratedBacteriaData.rds")
e1_test<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/Exp1","MH63","Root","ANCOMBC2_cutMin.rds",sep="_"))
tab_all1<-data.frame()
for(Exp in c("Exp1")){
	for(Com in c("Root","Rhizosphere")){
		for(Var in c("MH63","NIP")){
			e1<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin.rds",sep=""))
			if(length(which(e1$res[,"p_Timepointday3"] < 0.01))>0){
			e1_sig1<-e1$res[which(e1$res[,"p_Timepointday3"] < 0.01),]#$taxon
			e1_sig1_n<-e1_sig1[which(abs(e1_sig1[,"lfc_Timepointday3"]) > log2(1.5)),]#
			}else{
			e1_sig1_n<-c()
			}

			if(length(which(e1$res[,"p_Timepointday5"] < 0.01))>0){
			e1_sig2<-e1$res[which(e1$res[,"p_Timepointday5"] < 0.01),]#$taxon
			e1_sig2_n<-e1_sig2[which(abs(e1_sig2[,"lfc_Timepointday5"]) > log2(1.5)),]#
			}else{
			e1_sig2_n<-c()
			}

			if(length(which(e1$res[,"p_Timepointday7"] < 0.01))>0){
			e1_sig3<-e1$res[which(e1$res[,"p_Timepointday7"] < 0.01),]#$taxon
			e1_sig3_n<-e1_sig3[which(abs(e1_sig3[,"lfc_Timepointday7"]) > log2(1.5)),]#
			}else{
			e1_sig3_n<-c()
			}

			e1_sig_all<-rbind(e1_sig1_n,e1_sig2_n,e1_sig3_n)
			colnames(e1_sig_all)<-colnames(e1_test$res)
			tax<-data.frame(tax_table(IntData)[e1_sig_all$taxon,])
			tax$taxon<-rownames(tax)
			e1_sig_tax<-merge(e1_sig_all,tax, by.x = "taxon")
			e1_sig_tax$Experiment<-Exp
			e1_sig_tax$Compartment<-Com
			e1_sig_tax$Variety<-Var
			tab_all1<-rbind(tab_all1,e1_sig_tax)
			}
		}
	}

tab_all2<-data.frame()
for(Exp in c("Exp2")){
	for(Com in c("Root","Rhizosphere")){
		for(Var in c("MH63","NIP")){
			e1<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin.rds",sep=""))

			if(length(which(e1$res[,"p_Timepointday3"] < 0.01))>0){
			e1_sig1<-e1$res[which(e1$res[,"p_Timepointday3"] < 0.01),]#$taxon
			e1_sig1_n<-e1_sig1[which(abs(e1_sig1[,"lfc_Timepointday3"]) > log2(1.5)),]#
			}else{
			e1_sig1_n<-c()
			}
			if(length(which(e1$res[,"p_Timepointday4"] < 0.01))>0){
			e1_sig2<-e1$res[which(e1$res[,"p_Timepointday4"] < 0.01),]#$taxon
			e1_sig2_n<-e1_sig2[which(abs(e1_sig2[,"lfc_Timepointday4"]) > log2(1.5)),]#
			}else{
			e1_sig2_n<-c()
			}
			if(length(which(e1$res[,"p_Timepointday5"] < 0.01))>0){
			e1_sig3<-e1$res[which(e1$res[,"p_Timepointday5"] < 0.01),]#$taxon
			e1_sig3_n<-e1_sig3[which(abs(e1_sig3[,"lfc_Timepointday5"]) > log2(1.5)),]#
			}else{
			e1_sig3_n<-c()
			}
			e1_sig_all<-rbind(e1_sig1_n,e1_sig2_n,e1_sig3_n)
			colnames(e1_sig_all)<-colnames(e1_test$res)
			tax<-data.frame(tax_table(IntData)[e1_sig_all$taxon,])
			tax$taxon<-rownames(tax)
			e1_sig_tax<-merge(e1_sig_all,tax, by.x = "taxon")
			e1_sig_tax$Experiment<-Exp
			e1_sig_tax$Compartment<-Com
			e1_sig_tax$Variety<-Var
			tab_all2<-rbind(tab_all2,e1_sig_tax)
			}
		}
	}
tab_all<-rbind(tab_all1,tab_all2)
write.csv(tab_all,"OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_time_138.csv")




### significant result for defoliation x time group
rm(list = ls())
IntData = readRDS("IntegratedBacteriaData.rds")
tab_all1<-data.frame()
e1_test<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/","Exp1","_","MH63","_","Root","_","ANCOMBC2_cutMin_inter.rds",sep=""))
for(Exp in c("Exp1")){
	for(Com in c("Root","Rhizosphere")){
		for(Var in c("MH63","NIP")){
			e1<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin_inter.rds",sep=""))
			if(length(which(e1$res[,"p_Timepointday3:TreatmentIntact"] < 0.01))>0){
			e1_sig1<-e1$res[which(e1$res[,"p_Timepointday3:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig1_n<-e1_sig1[which(abs(e1_sig1[,"lfc_Timepointday3:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig1_n<-c()
			}

			if(length(which(e1$res[,"p_Timepointday5:TreatmentIntact"] < 0.01))>0){
			e1_sig2<-e1$res[which(e1$res[,"p_Timepointday5:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig2_n<-e1_sig2[which(abs(e1_sig2[,"lfc_Timepointday5:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig2_n<-c()

			}

			if(length(which(e1$res[,"p_Timepointday7:TreatmentIntact"] < 0.01))>0){
			e1_sig3<-e1$res[which(e1$res[,"p_Timepointday7:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig3_n<-e1_sig3[which(abs(e1_sig3[,"lfc_Timepointday7:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig3_n<-c()
			}
			e1_sig_all<-rbind(e1_sig1_n,e1_sig2_n,e1_sig3_n)

			colnames(e1_sig_all)<-colnames(e1_test$res)
			tax<-data.frame(tax_table(IntData)[e1_sig_all$taxon,])
			tax$taxon<-rownames(tax)

			e1_sig_tax<-merge(e1_sig_all,tax, by.x = "taxon")

			e1_sig_tax$Experiment<-Exp
			e1_sig_tax$Compartment<-Com
			e1_sig_tax$Variety<-Var
			tab_all1<-rbind(tab_all1,e1_sig_tax)
			}
		}
	}

tab_all2<-data.frame()
for(Exp in c("Exp2")){
	for(Com in c("Root","Rhizosphere")){
		for(Var in c("MH63","NIP")){
			e1<-readRDS(paste("OutputDir/ANCCOM-BC2_Figure4A/",Exp,"_",Var,"_",Com,"_","ANCOMBC2_cutMin_inter.rds",sep=""))
			if(length(which(e1$res[,"p_Timepointday3:TreatmentIntact"] < 0.01))>0){
			e1_sig1<-e1$res[which(e1$res[,"p_Timepointday3:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig1_n<-e1_sig1[which(abs(e1_sig1[,"lfc_Timepointday3:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig1_n<-c()

			}

			if(length(which(e1$res[,"p_Timepointday4:TreatmentIntact"] < 0.01))>0){
			e1_sig2<-e1$res[which(e1$res[,"p_Timepointday4:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig2_n<-e1_sig2[which(abs(e1_sig2[,"lfc_Timepointday4:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig2_n<-c()

			}

			if(length(which(e1$res[,"p_Timepointday5:TreatmentIntact"] < 0.01))>0){
			e1_sig3<-e1$res[which(e1$res[,"p_Timepointday5:TreatmentIntact"] < 0.01),]#$taxon
			e1_sig3_n<-e1_sig3[which(abs(e1_sig3[,"lfc_Timepointday5:TreatmentIntact"]) > log2(1.5)),]#
			}else{
			e1_sig3_n<-c()

			}
			e1_sig_all<-rbind(e1_sig1_n,e1_sig2_n,e1_sig3_n)
			colnames(e1_sig_all)<-colnames(e1_test$res)
			tax<-data.frame(tax_table(IntData)[e1_sig_all$taxon,])
			tax$taxon<-rownames(tax)
			e1_sig_tax<-merge(e1_sig_all,tax, by.x = "taxon")
			e1_sig_tax$Experiment<-Exp
			e1_sig_tax$Compartment<-Com
			e1_sig_tax$Variety<-Var
			tab_all2<-rbind(tab_all2,e1_sig_tax)
			}
		}
	}

tab_all<-rbind(tab_all1,tab_all2)


write.csv(tab_all,"OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_timexdefoliation_138.csv")


#### bar plot for all results


library(ggplot2)
tab0<-read.csv("OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_138.csv")
tab1<-read.csv("OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_timexdefoliation_138.csv")
tab2<-read.csv("OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_time_138.csv")

p<-list()
i<-1
for(Exp in c("Exp1","Exp2")){
    for(Var in c("NIP","MH63")){
		for(Com in c("Root","Rhizosphere")){

        diff_def1<-subset(tab0,Experiment == Exp)
        diff_def2<-subset(diff_def1,Compartment == Com)
        diff_def3<-subset(diff_def2,Variety == Var)
        
        diff_int1<-subset(tab1,Experiment == Exp)
        diff_int2<-subset(diff_int1,Compartment == Com)
        diff_int3<-subset(diff_int2,Variety == Var)
        
        
        diff_tim1<-subset(tab2,Experiment == Exp)
        diff_tim2<-subset(diff_tim1,Compartment == Com)
        diff_tim3<-subset(diff_tim2,Variety == Var)
        
        tab<-data.frame(treat=c("defoliation","defoliation x timepoint","timepoint"),count=c(length(unique(diff_def3$taxon)),length(unique(diff_int3$taxon)),length(unique(diff_tim3$taxon))))
        p[[i]]<-ggplot(tab, aes(x = treat, y = count,fill=treat)) + geom_col(width=0.5)+ylab("Number of functions")+#ylim(0,80)+
          theme(axis.text.x=element_blank(),axis.title.x=element_blank())+
          geom_text(   
            aes(label=count, 
                vjust = -0.5,
                hjust = 0.5, 
            ),
            size=3.5   
          )+
          scale_fill_manual(values =c("#E37933","#E37933","#406CB4"))+
          theme(legend.position = "none")
        
        i=i+1
        
      }
    }
  }
  
library(Rmisc)
  
pdf("OutputDir/ANCCOM-BC2_Figure4A/countZOTUs_3_factor.pdf",width=6,height=3)
multiplot(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],cols =4)
dev.off()
  
  

