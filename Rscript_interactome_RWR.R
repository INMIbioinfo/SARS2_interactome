################################################################################################################
################################################################################################################
####
#### Rscript_interactome_RWR: MODULE AND FUNCTION USED TO GENERATE THE HOST-VIRUS INTERACTOME.
#### 
################################################################################################################
################################################################################################################

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# many functions developed used functions of  Alberto Valdeolivas and colleagues (DOI:10.1093/bioinformatics/bty637) and reported in this file
# https://github.com/alberto-valdeolivas/RWR-MH/blob/master/Scripts_and_Files/Network_Generation/Network_Generation_Utils.R
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 1. RWR for single and multiple viral seeds in PPi network

setwd("https://github.com/INMIbioinfo/SARS2_interactome/")
source("./script/Network_Generation_Utils.R")
PPI_table <- read.table("./files/PPI_2020-02-25_07_13_SARS_CoV_2.gr")
library(igraph)
library(PSICQUIC)
library(biomaRt)
library(digest)
library(RandomWalkRestartMH)
library(graphite)
library(org.Hs.eg.db)
library(RandomWalkRestartMH)

SARS2 <- c("E","M","N","NS7b","nsp1","nsp10","nsp11","nsp12","nsp13","nsp14","nsp15","nsp16","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","ORF10","ORF14","ORF1a","ORF3a","ORF6","ORF7a","ORF8","ORF9b","S")
Seed = SARS2
Virus <- c("SARS-CoV-2")
dir.create("./SARS2_Interactome/")
setwd("./SARS2_Interactome/")

RWR_Virus<- function(X, PPI_table){
    	PPI_Network <- graph.data.frame(PPI_table,directed=FALSE)
    	PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
    	PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
    	PPI_MultiplexObject
    	AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
    	AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
    	for (SeedGene in X) {
	    RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_MultiplexObject,SeedGene)
	    RWR_PPI_Results
	    dir<-paste("./", SeedGene, "") 
	    dir.create(dir)
	    setwd(dir)
	    RWRM<-(RWR_PPI_Results$RWRM_Results)
	    prox<-paste("./", SeedGene, "_RWR_PPI_Results_proximity.txt") 
   	    write.table(RWRM, file=prox)
	    TopResults_PPI <-create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,k=50)
 	    par(mar=c(0.1,0.1,0.1,0.1))
	    network_raw<-paste("./", SeedGene, "_RWR_PPI_new_Results_network_raw.tif") 
 	    tiff(file=network_raw)
 	    plot(TopResults_PPI, vertex.label.color="black", vertex.label.cex = .75, vertex.frame.color="#ffffff", vertex.size= 5, edge.curved=.2, vertex.color = ifelse(igraph::V(TopResults_PPI)$name == SeedGene,"yellow","#00CCFF"), edge.color="blue",edge.width=0.8, mark.groups = SeedGene, mark.col = "red", mark.border = NA, mark.expand = 15)
 	    dev.off()
 	    top<-as_data_frame(TopResults_PPI)
 	    vert_tor<-as_data_frame(TopResults_PPI, what= "vertices")
	    network<-paste("./", SeedGene, "_RWR_PPI_new_Results_network.txt") 
 	    write.table(top, file=network)
	    ver<-paste("./", SeedGene, "_RWR_PPI_new_Results_ver.txt") 
 	    write.table(vert_tor, file=ver)
	    setwd("../.")
	    }
	RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_MultiplexObject,Seed)
	RWR_PPI_Results
	dir<-paste("./", Virus, "") 
	dir.create(dir)
	setwd(dir)
	RWRM<-(RWR_PPI_Results$RWRM_Results)
	prox<-paste("./", Virus, "_RWR_PPI_Results_proximity.txt") 
   	write.table(RWRM, file=prox)
	TopResults_PPI <-create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,k=199)
 	par(mar=c(0.1,0.1,0.1,0.1))
	network_raw<-paste("./", Virus, "_RWR_PPI_new_Results_network_raw.tif") 
 	tiff(file=network_raw)
 	plot(TopResults_PPI, vertex.label.color="black", vertex.label.cex = .75, vertex.frame.color="#ffffff", vertex.size= 5, edge.curved=.2, vertex.color = ifelse(igraph::V(TopResults_PPI)$name == "S","yellow","#00CCFF"), edge.color="blue",edge.width=0.8, mark.groups = Seed, mark.col = "red", mark.border = NA, mark.expand = 15)
 	dev.off()
 	top<-as_data_frame(TopResults_PPI)
 	vert_tor<-as_data_frame(TopResults_PPI, what= "vertices")
	network<-paste("./", Virus, "_RWR_PPI_new_Results_network.txt") 
 	write.table(top, file=network)
	ver<-paste("./", Virus, "_RWR_PPI_new_Results_ver.txt") 
 	write.table(vert_tor, file=ver)
	setwd("../.")
	  }
RWR_Virus(Seed, PPI_table)





