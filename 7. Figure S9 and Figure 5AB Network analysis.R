
rm(list = ls())

###1. intact-netwrok & defoliation-network


library(phyloseq)
library(Hmisc)
IntData = readRDS("IntegratedBacteriaData.rds")     

for(Tre in c("Intact","Defoliation")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			IntDataSelTre =subset_samples(IntData, Treatment %in% Tre)
			IntDataSelTreVar =subset_samples(IntDataSelTre, Variety %in% Var)
			IntDataSelTreVarCom<- subset_samples(IntDataSelTreVar, Compartment %in% Com)
			
			IntDataSelTreVarComNor <- transform_sample_counts(IntDataSelTreVarCom, function(x) x / sum(x) )
			IntDataSelTreVarComNorFilt = filter_taxa(IntDataSelTreVarComNor, function(x) mean(x) > 0.0001, TRUE)
			IntDataSelTreVarComNorFilt = filter_taxa(IntDataSelTreVarComNorFilt,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
			
			meta<-sample_data(IntDataSelTreVarCom)
			otu_table<-data.frame(otu_table(IntDataSelTreVarCom))
			count =  data.frame(otu_table)
			library(edgeR)
			Treat<- factor(meta$Treatment)
			Time<- factor(substring(meta$Timepoint,4,4))
			y <- DGEList(counts=count, group=Treat)
			y <- calcNormFactors(y)
			zotutab<- cpm(y)
			zotutab<-zotutab[rownames(otu_table(IntDataSelTreVarComNorFilt)),]
			comm.data<-data.frame(t(data.frame(zotutab)))
			occor<-rcorr(as.matrix(comm.data),type="spearman")
			r.cor = occor$r 
			p.cor = occor$P 
			p.adj <- p.adjust(p.cor, method="BH")
			r.matrix<-r.cor
			p.matrix<-p.cor
			m1<-r.matrix
			m2<-p.matrix
			for(i in 1:nrow(m1)){
				m1[i,i:ncol(m1)]<-rep(0,length( m1[i,i:ncol(m1)]))
				}
			network_m1<-data.frame()
			for(i in 1:nrow(m1)){
				if(length(colnames(m1)[which(abs(m1[i,])>0)])>0){
					tab<-data.frame(taxa1=rep(rownames(m1)[i],length(colnames(m1)[which(abs(m1[i,])>0)])),taxa2=colnames(m1)[which(abs(m1[i,])>0)],r=m1[i,][which(abs(m1[i,])>0)],p=m2[i,][which(abs(m1[i,])>0)])
					rownames(tab)<-paste(tab$taxa1,tab$taxa2,sep="_")
					network_m1<-rbind(network_m1,tab)
					}
				}
        

			#network_m1<-tab0_bs
			colnames(network_m1)[1:2]<-c("Source","Target")
			network_m1<-network_m1[which(network_m1$r != 1),]
			network_m1$adjust.p<-p.adjust(network_m1$p, method="BH")
			network_m1<-subset(network_m1,abs(r)>0.8)
			network_m1<-subset(network_m1,adjust.p<0.001)
			name<-unique(c(as.character(network_m1$Source),as.character(network_m1$Target)))
			network_m2<-data.frame(id=name,label=name)
			tax<-data.frame(tax_table(IntDataSelTreVarComNorFilt))
			tax$label<-rownames(tax)
			network_m2<- merge(network_m2,tax, by.x = "label")
			otu<-data.frame(otu_table(IntDataSelTreVarComNorFilt))
			table<-data.frame(name=rownames(otu),abun=rowMeans(otu))
			table$label<-rownames(table)
			network_m2<- merge(network_m2,table, by.x = "label")
			write.csv(network_m1,paste("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/net",Tre,Var,Com,"network_edges",".csv",sep="_"),row.names=FALSE)
			write.csv(network_m2,paste("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/net",Tre,Var,Com,"network_nodes",".csv",sep="_"),row.names=FALSE)
        
			}
		}
	}





###The difference in degree between the int and def networks (degree calculated using Gephi).



library(ggplot2)	
library(ggsignif)
tab_all<-data.frame()
for(Tre in c("Intact","Defoliation")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			tab<-read.csv(paste("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/addDegree/net",Tre,Var,Com,"network_nodes_AddDegree.csv",sep="_"))
			tab$GroupMergeVarCom=paste(Var,Com,"NODEs",sep="_")
			tab$GroupTre=Tre
			tab_all<-rbind(tab_all,tab)
			}
		}
	}

tab_all$GroupTre<-factor(tab_all$GroupTre,level=c("Intact","Defoliation"))
tab_all$GroupMergeVarCom<-factor(tab_all$GroupMergeVarCom,level=c("MH63_Root_NODEs","MH63_Rhizosphere_NODEs","NIP_Root_NODEs","NIP_Rhizosphere_NODEs"))
p<-ggplot(tab_all, aes(x=GroupTre,y=degree,group=GroupTre,fill=GroupTre)) +
	geom_jitter(aes(fill=GroupTre,color=GroupTre),position = position_jitter(width = 0.2))+facet_grid(.~GroupMergeVarCom)+
    geom_signif(
    comparisons = list(
     c("Intact","Defoliation")
    ), #检测两者之间的差异显著性
    map_signif_level = T, #添加星号标记
    test = "wilcox.test", #检测方法
    vjust=0.1, #标注和横线的距离
    tip_length = 0.05 #两端短竖线的长度
	)+
	scale_colour_manual(values = c("#406CB4","#E37933"))

pdf("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/all_degree_comparision12x5.pdf",width=12,height=5)
p
dev.off()

###The difference in the number of nodes between the int and def networks



tab_all<-c()
for(Tre in c("Intact","Defoliation")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			tab<-read.csv(paste("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/addDegree/net",Tre,Var,Com,"network_nodes_AddDegree.csv",sep="_"))
			t1<-c(nrow(tab),Tre,paste(Var,Com,"NODEs",sep="_"))
			tab_all<-rbind(tab_all,t1)
			}
		}
	}
	
	
	
	
colnames(tab_all)<-c("count","GroupTre","GroupMergeVarCom")
tab_all<-data.frame(tab_all)
tab_all$count<-as.numeric(tab_all$count)
tab_all$GroupTre<-factor(tab_all$GroupTre,level=c("Intact","Defoliation"))
tab_all$GroupMergeVarCom<-factor(tab_all$GroupMergeVarCom,level=c("MH63_Root_NODEs","MH63_Rhizosphere_NODEs","NIP_Root_NODEs","NIP_Rhizosphere_NODEs"))

p<-ggplot(tab_all, aes(x=GroupTre, y=count, fill=GroupTre)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.5)+facet_grid(.~GroupMergeVarCom)+
		    geom_text(aes(label=count),size=4,vjust=-0.5)+ylim(0,1200)+
			 scale_fill_manual(values = c("#406CB4","#E37933"))

pdf("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/all_nodes_comparision12x5_width5.pdf",width=12,height=5)
p
dev.off()


###The difference in the number of edges between the int and def networks


tab_all<-c()
for(Tre in c("Intact","Defoliation")){
	for(Var in c("MH63","NIP")){
		for(Com in c("Root","Rhizosphere")){
			tab<-read.csv(paste("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/net",Tre,Var,Com,"network_edge_m",".csv",sep="_"))
			t1<-c(nrow(tab),Tre,paste(Var,Com,"NODEs",sep="_"))
			tab_all<-rbind(tab_all,t1)
			}
		}
	}
colnames(tab_all)<-c("count","GroupTre","GroupMergeVarCom")
tab_all<-data.frame(tab_all)
tab_all$count<-as.numeric(tab_all$count)
tab_all$GroupTre<-factor(tab_all$GroupTre,level=c("Intact","Defoliation"))
tab_all$GroupMergeVarCom<-factor(tab_all$GroupMergeVarCom,level=c("MH63_Root_NODEs","MH63_Rhizosphere_NODEs","NIP_Root_NODEs","NIP_Rhizosphere_NODEs"))
p<-ggplot(tab_all, aes(x=GroupTre, y=count, fill=GroupTre)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.5)+facet_grid(.~GroupMergeVarCom)+
		    geom_text(aes(label=count),size=4,vjust=-0.5)+ylim(0,27000)+
			 scale_fill_manual(values = c("#406CB4","#E37933"))



pdf("OutputDir/NetworkAnalysis/IntactAndDefoliationNetwork/all_edges_comparision12x5_width5.pdf",width=12,height=5)
p
dev.off()






###2. Root-network & Rhizo-network

##Root
library(edgeR)
library(phyloseq)
library(Hmisc)
DiffZotuTable = read.csv("OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_138.csv")
DiffZotuTableRe = subset(DiffZotuTable,Compartment == "Root")
IntData = readRDS("IntegratedBacteriaData.rds")
IntData = subset_samples(IntData, Treatment %in% c("Intact"))
IntData = subset_samples(IntData, Compartment %in% c("Root"))
SampleData<-sample_data(IntData)
OtuTable<-data.frame(otu_table(IntData))
Treat<- factor(SampleData$Treatment)
Time<- factor(substring(SampleData$Timepoint,4,4))
y <- DGEList(counts=OtuTable, group=Treat)
y <- calcNormFactors(y)
OtuTableCpm<- cpm(y)
IntDataSel = transform_sample_counts(IntData, function(x) x / sum(x) )
IntDataSel = filter_taxa(IntDataSel, function(x) mean(x) > 0.0001, TRUE)
IntDataSel = filter_taxa(IntDataSel,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
OtuTableCpm<-OtuTableCpm[rownames(otu_table(IntDataSel)),]
NetworkData<-data.frame(t(data.frame(OtuTableCpm)))
occor<-rcorr(as.matrix(NetworkData),type="spearman")
r.cor = occor$r 
p.cor = occor$P 
p.adj <- p.adjust(p.cor, method="BH")
r.matrix<-r.cor
p.matrix<-p.cor
for(i in 1:nrow(r.matrix)){
    r.matrix[i,i:ncol(r.matrix)]<-rep(0,length( r.matrix[i,i:ncol(r.matrix)]))
    }
EdgeTable<-data.frame()
for(i in 1:nrow(r.matrix)){
    if(length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])>0){
       tab<-data.frame(Source=rep(rownames(r.matrix)[i],length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])),Target=colnames(r.matrix)[which(abs(r.matrix[i,])>0)],r=r.matrix[i,][which(abs(r.matrix[i,])>0)],p=p.matrix[i,][which(abs(r.matrix[i,])>0)])
       rownames(tab)<-paste(tab$Source,tab$Target,sep="_")
       EdgeTable<-rbind(EdgeTable,tab)
       }
    }
EdgeTable<-EdgeTable[which(EdgeTable$r != 1),]
EdgeTable$adjust.p<-p.adjust(EdgeTable$p, method="BH")
EdgeTable<-subset(EdgeTable,abs(r)>0.8)
EdgeTable<-subset(EdgeTable,adjust.p<0.001)
name<-unique(c(as.character(EdgeTable$Source),as.character(EdgeTable$Target)))
NodesTable<-data.frame(id=name,label=name)
NodesTable$labelDiff<-"other"
NodesTable[which(NodesTable$label %in% DiffZotuTableRe$taxon),]$labelDiff<-rep("diff",length(NodesTable[which(NodesTable$label %in% DiffZotuTableRe$taxon),]$labelDiff))
TaxTable<-data.frame(tax_table(IntDataSel))
TaxTable$label<-rownames(TaxTable)
NodesTableAddTax<- merge(NodesTable,TaxTable, by.x = "label")
OtuTableSel<-data.frame(otu_table(IntDataSel))
OtuAbunMatrix<-data.frame(name=rownames(OtuTableSel),abun=rowMeans(OtuTableSel))
OtuAbunMatrix$label<-rownames(OtuAbunMatrix)
NodesTableAddTaxAbun<- merge(NodesTableAddTax,OtuAbunMatrix, by.x = "label")
write.csv(EdgeTable,paste("OutputDir/NetworkAnalysis/Edges_Table_Root_network",".csv",sep=""),row.names=FALSE)
write.csv(NodesTableAddTaxAbun,paste("OutputDir/NetworkAnalysis/Nodes_Table_Root_network",".csv",sep=""),row.names=FALSE)





        
##Rhizosphere
library(edgeR)
library(phyloseq)
library(Hmisc)
DiffZotuTable = read.csv("OutputDir/ANCCOM-BC2_Figure4A/ANCOMBC_def_result_138.csv")
DiffZotuTableRs = subset(DiffZotuTable,Compartment == "Rhizosphere")
IntData = readRDS("IntegratedBacteriaData.rds")
IntData = subset_samples(IntData, Treatment %in% c("Intact"))
IntData = subset_samples(IntData, Compartment %in% c("Rhizosphere"))
SampleData<-sample_data(IntData)
OtuTable<-data.frame(otu_table(IntData))
Treat<- factor(SampleData$Treatment)
Time<- factor(substring(SampleData$Timepoint,4,4))
y <- DGEList(counts=OtuTable, group=Treat)
y <- calcNormFactors(y)
OtuTableCpm<- cpm(y)
IntDataSel = transform_sample_counts(IntData, function(x) x / sum(x) )
IntDataSel = filter_taxa(IntDataSel, function(x) mean(x) > 0.0001, TRUE)
IntDataSel = filter_taxa(IntDataSel,function(x) sum(x > 0) > (0.1*length(x)),TRUE)
OtuTableCpm<-OtuTableCpm[rownames(otu_table(IntDataSel)),]
NetworkData<-data.frame(t(data.frame(OtuTableCpm)))
occor<-rcorr(as.matrix(NetworkData),type="spearman")
r.cor = occor$r # 取相关性矩阵R值
p.cor = occor$P # 取相关性矩阵p值
p.adj <- p.adjust(p.cor, method="BH")
r.matrix<-r.cor
p.matrix<-p.cor
for(i in 1:nrow(r.matrix)){
r.matrix[i,i:ncol(r.matrix)]<-rep(0,length( r.matrix[i,i:ncol(r.matrix)]))
}
EdgeTable<-data.frame()
for(i in 1:nrow(r.matrix)){
    if(length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])>0){
       tab<-data.frame(Source=rep(rownames(r.matrix)[i],length(colnames(r.matrix)[which(abs(r.matrix[i,])>0)])),Target=colnames(r.matrix)[which(abs(r.matrix[i,])>0)],r=r.matrix[i,][which(abs(r.matrix[i,])>0)],p=p.matrix[i,][which(abs(r.matrix[i,])>0)])
       rownames(tab)<-paste(tab$Source,tab$Target,sep="_")
       EdgeTable<-rbind(EdgeTable,tab)
       }
    }	
EdgeTable<-EdgeTable[which(EdgeTable$r != 1),]
EdgeTable$adjust.p<-p.adjust(EdgeTable$p, method="BH")
EdgeTable<-subset(EdgeTable,abs(r)>0.8)
EdgeTable<-subset(EdgeTable,adjust.p<0.001)
name<-unique(c(as.character(EdgeTable$Source),as.character(EdgeTable$Target)))
NodesTable<-data.frame(id=name,label=name)
NodesTable$labelDiff<-"other"
NodesTable[which(NodesTable$label %in% DiffZotuTableRs$Row.names),]$labelDiff<-rep("diff",length(NodesTable[which(NodesTable$label %in% DiffZotuTableRs$Row.names),]$labelDiff))
TaxTable<-data.frame(tax_table(IntDataSel))
TaxTable$label<-rownames(TaxTable)
NodesTableAddTax<- merge(NodesTable,TaxTable, by.x = "label")
OtuTableSel<-data.frame(otu_table(IntDataSel))
OtuAbunMatrix<-data.frame(name=rownames(OtuTableSel),abun=rowMeans(OtuTableSel))
OtuAbunMatrix$label<-rownames(OtuAbunMatrix)
NodesTableAddTaxAbun<- merge(NodesTableAddTax,OtuAbunMatrix, by.x = "label")
write.csv(EdgeTable,paste("OutputDir/NetworkAnalysis/Edges_Table_Rhizo_network",".csv",sep=""),row.names=FALSE)
write.csv(NodesTableAddTaxAbun,paste("OutputDir/NetworkAnalysis/Nodes_Table_Rhizo_network",".csv",sep=""),row.names=FALSE)



























		