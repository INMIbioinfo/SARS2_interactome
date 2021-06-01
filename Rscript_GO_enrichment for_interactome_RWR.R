# 2. Gene ontology enrichment on Enrichr website

library(enrichR)

GO_Virus<- function(X){
	for (SeedGene in X) {
	    dir<-paste("./", SeedGene, "") 
 	    setwd(dir)
 	    dbs_1 <- c("WikiPathways_2019_Human")
 	    dbs_2 <- c("COVID-19_Related_Gene_Sets")
 	    dbs_3 <- c("Reactome_2016")
 	    dbs_4 <- c("KEGG_2019_Human")
	    ver<-paste("./", SeedGene, "_RWR_PPI_new_Results_ver.txt") 
 	    vert_tor <- read.table(file=ver)
 	    colnames(vert_tor)<- NULL
 	    rownames(vert_tor)<- NULL
 	    vert_tor[,2]<- NULL
 	    vert_tor<-as.data.frame(vert_tor)
 	    vert_tor<-as.character(vert_tor[,1])
 	    enriched_1 <- enrichr(vert_tor, dbs_1)
 	    enriched_2 <- enrichr(vert_tor, dbs_2)
 	    enriched_3 <- enrichr(vert_tor, dbs_3)
 	    enriched_4 <- enrichr(vert_tor, dbs_4)
 	    wiki<-paste("./", SeedGene, "_RWR_PPI_new_Results_GO_enrichment_WikiPathways_2019_Human.txt")
 	    covid<-paste("./", SeedGene, "_RWR_PPI_new_Results_GO_enrichment_COVID-19_Related_Gene_Sets.txt")
 	    reactome<-paste("./", SeedGene, "_RWR_PPI_new_Results_GO_enrichment_Reactome_2016.txt")
 	    kegg<-paste("./", SeedGene, "_RWR_PPI_new_Results_GO_enrichment_KEGG_2019_Human.txt")
 	    write.table(enriched_1, file=wiki)
 	    write.table(enriched_2, file=covid)
 	    write.table(enriched_3, file=reactome)
 	    write.table(enriched_4, file=kegg)
	    setwd("../.")
	    }
}
GO_Virus(Seed)
