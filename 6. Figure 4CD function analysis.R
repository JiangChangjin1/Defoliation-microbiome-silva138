



#####1.differential functions analysis

library(ggpicrust2)
library(phyloseq)

# Convert the PICRUSt2-based predicted KO table to a KEGG pathway table.
ko_abundance_file <- "OutputDir/Picrust2_predicted_kegg/pred_metagenome_unstrat.tsv"
kegg_abundanceAbun <- ko2kegg_abundance(file = ko_abundance_file)


IntData  = readRDS("IntegratedBacteriaData.rds")
AllExp1Rel<-c()
for(Com in c("Root","Rhizosphere")){
	for(Var in c("NIP","MH63")){
		for(Exp in c("Exp1")){
			for(Tim in c("day1","day3","day5","day7")){
				IntDataSelectCom<- subset_samples(IntData , Compartment %in% Com)
				IntDataSelectComExp<- subset_samples(IntDataSelectCom, Experiment %in% Exp)
				IntDataSelectComExpVar<- subset_samples(IntDataSelectComExp, Variety %in% Var )
				IntDataSelectComExpVarTim<- subset_samples(IntDataSelectComExpVar, Timepoint %in% Tim )
				metadata <- data.frame(sample_data(IntDataSelectComExpVarTim))
				kegg_abundanceSel<-kegg_abundanceAbun[,metadata$SampleID ]
				if(length(which(rowSums(kegg_abundanceSel) == 0 ))>0){
					kegg_abundanceSel<-kegg_abundanceSel[-which(rowSums(kegg_abundanceSel) == 0 ),]
					}
				log2FC_all<-c()
				def_all<-c()
				int_all<-c()
				
				###normalized kegg_abundance table for calculating log2FC
				df<-t(kegg_abundanceSel)
				data.prop<-df/rowSums(df)
				kegg_abundanceSelNor<-t(data.prop)
				for(i in 1:nrow(kegg_abundanceSelNor)){
					log2FC<-log2(mean(t(kegg_abundanceSelNor)[,i][1:3])/mean(t(kegg_abundanceSelNor)[,i][4:6]))
					log2FC_all<-c(log2FC_all,log2FC)
					def<-mean(t(kegg_abundanceSelNor)[,i][1:3])
					int<-mean(t(kegg_abundanceSelNor)[,i][4:6])	
					def_all<-c(def_all,def)
					int_all<-c(int_all,int)
					}

				daa_results_df <- pathway_daa(abundance = kegg_abundanceSel, metadata = metadata, group = "Treatment", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = NULL)
				daa_results_df$Compartment <- Com
				daa_results_df$Experiment <- Exp
				daa_results_df$Variety <- Var
				daa_results_df$Timepoint <- Tim
				daa_results_df$logFC<-log2FC_all
				daa_results_df$int_all<-int_all
				daa_results_df$def_all<-def_all
				AllExp1Rel<-rbind(AllExp1Rel,daa_results_df)

				}


			}

		}
}





IntData  = readRDS("IntegratedBacteriaData.rds")
AllExp2Rel<-c()
for(Com in c("Root","Rhizosphere")){
	for(Var in c("NIP","MH63")){
		for(Exp in c("Exp2")){
			for(Tim in c("day2","day3","day4","day5")){
				IntDataSelectCom<- subset_samples(IntData , Compartment %in% Com)
				IntDataSelectComExp<- subset_samples(IntDataSelectCom, Experiment %in% Exp)
				IntDataSelectComExpVar<- subset_samples(IntDataSelectComExp, Variety %in% Var )
				IntDataSelectComExpVarTim<- subset_samples(IntDataSelectComExpVar, Timepoint %in% Tim )
				metadata <- data.frame(sample_data(IntDataSelectComExpVarTim))
				kegg_abundanceSel<-kegg_abundanceAbun[,metadata$SampleID ]
				if(length(which(rowSums(kegg_abundanceSel) == 0 ))>0){
					kegg_abundanceSel<-kegg_abundanceSel[-which(rowSums(kegg_abundanceSel) == 0 ),]
					}
				log2FC_all<-c()
				def_all<-c()
				int_all<-c()
				
				###normalized kegg_abundance table for calculating log2FC
				df<-t(kegg_abundanceSel)
				data.prop<-df/rowSums(df)
				kegg_abundanceSelNor<-t(data.prop)
				for(i in 1:nrow(kegg_abundanceSelNor)){
					log2FC<-log2(mean(t(kegg_abundanceSelNor)[,i][1:3])/mean(t(kegg_abundanceSelNor)[,i][4:6]))
					log2FC_all<-c(log2FC_all,log2FC)
					def<-mean(t(kegg_abundanceSelNor)[,i][1:3])
					int<-mean(t(kegg_abundanceSelNor)[,i][4:6])	
					def_all<-c(def_all,def)
					int_all<-c(int_all,int)
					}

				daa_results_df <- pathway_daa(abundance = kegg_abundanceSel, metadata = metadata, group = "Treatment", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = NULL)
				daa_results_df$Compartment <- Com
				daa_results_df$Experiment <- Exp
				daa_results_df$Variety <- Var
				daa_results_df$Timepoint <- Tim
				daa_results_df$logFC<-log2FC_all
				daa_results_df$int_all<-int_all
				daa_results_df$def_all<-def_all
				AllExp2Rel<-rbind(AllExp2Rel,daa_results_df)

				}


			}

		}
}




CombTwoExpRel<-rbind(AllExp1Rel,AllExp2Rel)
write.csv(CombTwoExpRel,"OutputDir/Picrust2_predicted_kegg/Combinded_TwoExp_ggpicrust2_LinDA_result.csv")





#####2.Add different kegg level infomation



  
rm(list = ls())
	  
	  
CombTwoExpRel<-read.csv("OutputDir/Picrust2_predicted_kegg/Combinded_TwoExp_ggpicrust2_LinDA_result.csv")
KeggPathway<-read.csv("OutputDir/Picrust2_predicted_kegg/KEGG-pathway-classification-m.csv")
KeggPathway<-KeggPathway[,c("feature","level1","level2","level3")]
CombTwoExpRel$level1<-"Unassigned"
CombTwoExpRel$level2<-"Unassigned"
CombTwoExpRel$level3<-"Unassigned"
for(i in 1:nrow(CombTwoExpRel)){
	if(length(which(KeggPathway$feature %in% CombTwoExpRel$feature[i]))>0){
		CombTwoExpRel$level1[i]<-KeggPathway$level1[which(KeggPathway$feature %in% CombTwoExpRel$feature[i])]
		CombTwoExpRel$level2[i]<-KeggPathway$level2[which(KeggPathway$feature %in% CombTwoExpRel$feature[i])]
		CombTwoExpRel$level3[i]<-KeggPathway$level3[which(KeggPathway$feature %in% CombTwoExpRel$feature[i])]
		}
	}
CombTwoExpRel<-CombTwoExpRel[which(CombTwoExpRel$level3 %in%  c("Metabolism","Environmental Information Processing")),]
write.csv(CombTwoExpRel,"OutputDir/Picrust2_predicted_kegg/all_function_ggpicrust2_ori_Abun_LinDA_add_level2Info.csv")
	  
	  



#####3. Screening for predicted functions that respond to defoliation with the same trend over two experiments.



rm(list = ls())

allKegg<-read.csv("OutputDir/Picrust2_predicted_kegg/all_function_ggpicrust2_ori_Abun_LinDA_add_level2Info.csv")
allKeggSig<-subset(allKegg,p_values<0.01)

allKeggSig$group<-paste(allKeggSig$level1 ,allKeggSig$Compartment ,allKeggSig$Experiment,allKeggSig$Variety,sep="_")
allKeggSigEnri<-allKeggSig[which(allKeggSig$logFC > 0),]
allKeggSigDep<-allKeggSig[which(allKeggSig$logFC < 0),]

###Remove results that are both enriched and depleted in the same group of time series.
delectEnri<-unique(allKeggSigEnri$group)[which(unique(allKeggSigEnri$group) %in% unique(allKeggSigDep$group))]
delectDep<-unique(allKeggSigDep$group)[which(unique(allKeggSigDep$group) %in% unique(allKeggSigEnri$group))]
allKeggSigEnri<-allKeggSigEnri[-which(allKeggSigEnri$group %in% delectEnri),]
allKeggSigDep<-allKeggSigDep[-which(allKeggSigDep$group %in% delectDep),]
write.csv(data.frame(rbind(allKeggSigEnri,allKeggSigDep)),"OutputDir/Picrust2_predicted_kegg/KeggSigEnrichedOrDepletedInDefoliation.csv")

#length(unique(c(allKeggSigEnri$feature,allKeggSigDep$feature)))
#83

### Identification of functions stably enriched over two experiments.

rm(list = ls())

allKeggSig<-allKeggSigEnri

allKeggSigNip<-subset(allKeggSig,Variety == "NIP")
allKeggSigNipRs<-subset(allKeggSigNip,Compartment  == "Rhizosphere")
allKeggSigNipRe<-subset(allKeggSigNip,Compartment  == "Root")

NipSigNameRsExp1<-unique(subset(allKeggSigNipRs,Experiment  == "Exp1")$feature)
NipSigNameRsExp2<-unique(subset(allKeggSigNipRs,Experiment  == "Exp2")$feature) 
NipSigNameRsExp1[which(NipSigNameRsExp1 %in% NipSigNameRsExp2)]
character(0)

NipSigNameReExp1<-unique(subset(allKeggSigNipRe,Experiment  == "Exp1")$feature)
NipSigNameReExp2<-unique(subset(allKeggSigNipRe,Experiment  == "Exp2")$feature)
NipSigNameReExp1[which(NipSigNameReExp1 %in% NipSigNameReExp2)]
character(0)

allKeggSigMH63<-subset(allKeggSig,Variety == "MH63")
allKeggSigMH63Rs<-subset(allKeggSigMH63,Compartment  == "Rhizosphere")
allKeggSigMH63Re<-subset(allKeggSigMH63,Compartment  == "Root")

MH63SigNameRsExp1<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp1")$feature)
MH63SigNameRsExp2<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp2")$feature)
MH63SigNameRsExp1[which(MH63SigNameRsExp1 %in% MH63SigNameRsExp2)]
character(0)

MH63SigNameReExp1<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp1")$feature)
MH63SigNameReExp2<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp2")$feature)
MH63SigNameReExp1[which(MH63SigNameReExp1 %in% MH63SigNameReExp2)]
character(0)


length(NipSigNameReExp1)

length(NipSigNameReExp2)


length(NipSigNameRsExp1)

length(NipSigNameRsExp2)

length(MH63SigNameReExp1)

length(MH63SigNameReExp2)


length(MH63SigNameRsExp1)

length(MH63SigNameRsExp2)



### Identification of functions stably depleted over two experiments.

rm(list = ls())

allKeggSig<-allKeggSigDep

allKeggSigNip<-subset(allKeggSig,Variety == "NIP")
allKeggSigNipRs<-subset(allKeggSigNip,Compartment  == "Rhizosphere")
allKeggSigNipRe<-subset(allKeggSigNip,Compartment  == "Root")

NipSigNameRsExp1<-unique(subset(allKeggSigNipRs,Experiment  == "Exp1")$feature)
NipSigNameRsExp2<-unique(subset(allKeggSigNipRs,Experiment  == "Exp2")$feature)
NipSigNameRsExp1[which(NipSigNameRsExp1 %in% NipSigNameRsExp2)]
character(0)

NipSigNameReExp1<-unique(subset(allKeggSigNipRe,Experiment  == "Exp1")$feature)
NipSigNameReExp2<-unique(subset(allKeggSigNipRe,Experiment  == "Exp2")$feature)
NipSigNameReExp1[which(NipSigNameReExp1 %in% NipSigNameReExp2)]


allKeggSigMH63<-subset(allKeggSig,Variety == "MH63")
allKeggSigMH63Rs<-subset(allKeggSigMH63,Compartment  == "Rhizosphere")
allKeggSigMH63Re<-subset(allKeggSigMH63,Compartment  == "Root")
character(0)

MH63SigNameRsExp1<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp1")$feature)
MH63SigNameRsExp2<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp2")$feature)
MH63SigNameRsExp1[which(MH63SigNameRsExp1 %in% MH63SigNameRsExp2)]

MH63SigNameReExp1<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp1")$feature)
MH63SigNameReExp2<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp2")$feature)
MH63SigNameReExp1[which(MH63SigNameReExp1 %in% MH63SigNameReExp2)]

> MH63SigNameRsExp1[which(MH63SigNameRsExp1 %in% MH63SigNameRsExp2)]
[1] "ko00590" "ko00625" "ko00627" "ko00626" "ko00830" "ko00903" "ko00930"


> MH63SigNameReExp1[which(MH63SigNameReExp1 %in% MH63SigNameReExp2)]
 [1] "ko00310" "ko00380" "ko00627" "ko00626" "ko00071" "ko00072" "ko00363" "ko02010" "ko00040" "ko01040"
[11] "ko00280" "ko00640" "ko00410" "ko00650" "ko00903" "ko00930"

length(NipSigNameReExp1)

length(NipSigNameReExp2)


length(NipSigNameRsExp1)

length(NipSigNameRsExp2)

length(MH63SigNameReExp1)

length(MH63SigNameReExp2)


length(MH63SigNameRsExp1)

length(MH63SigNameRsExp2)






allKeggSigMH63<-subset(allKeggSig,Variety == "MH63")
allKeggSigMH63Rs<-subset(allKeggSigMH63,Compartment  == "Rhizosphere")
allKeggSigMH63Re<-subset(allKeggSigMH63,Compartment  == "Root")


MH63SigNameRsExp1_l<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp1")$level1)
MH63SigNameRsExp2_l<-unique(subset(allKeggSigMH63Rs,Experiment  == "Exp2")$level1)
MH63SigNameRsExp1_l[which(MH63SigNameRsExp1_l %in% MH63SigNameRsExp2_l)]

1] "Arachidonic acid metabolism"              
[2] "Chloroalkane and chloroalkene degradation"
[3] "Aminobenzoate degradation"                
[4] "Naphthalene degradation"                  
[5] "Retinol metabolism"                       
[6] "Limonene and pinene degradation"          
[7] "Caprolactam degradation"  

MH63SigNameReExp1_l<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp1")$level1)
MH63SigNameReExp2_l<-unique(subset(allKeggSigMH63Re,Experiment  == "Exp2")$level1)
MH63SigNameReExp1_l[which(MH63SigNameReExp1_l %in% MH63SigNameReExp2_l)]

 [1] "Lysine degradation"                        
 [2] "Tryptophan metabolism"                     
 [3] "Aminobenzoate degradation"                 
 [4] "Naphthalene degradation"                   
 [5] "Fatty acid degradation"                    
 [6] "Synthesis and degradation of ketone bodies"
 [7] "Bisphenol degradation"                     
 [8] "ABC transporters"                          
 [9] "Pentose and glucuronate interconversions"  
[10] "Biosynthesis of unsaturated fatty acids"   
[11] "Valine, leucine and isoleucine degradation"
[12] "Propanoate metabolism"                     
[13] "beta-Alanine metabolism"                   
[14] "Butanoate metabolism"                      
[15] "Limonene and pinene degradation"           
[16] "Caprolactam degradation" 

MH63CommRe_l<-MH63SigNameReExp1_l[which(MH63SigNameReExp1_l %in% MH63SigNameReExp2_l)]
MH63CommRs_l<-MH63SigNameRsExp1_l[which(MH63SigNameRsExp1_l %in% MH63SigNameRsExp2_l)]

saveRDS(MH63CommRe_l,"OutputDir/Picrust2_predicted_kegg/Def_MH63CommRe_depleted.rds")

saveRDS(MH63CommRs_l,"OutputDir/Picrust2_predicted_kegg/Def_MH63CommRs_depleted.rds")





