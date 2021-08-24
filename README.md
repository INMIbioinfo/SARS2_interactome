# SARS2_interactome

a network analysis on protein–protein interactions (PPI) between viral and host proteins to better identify host biological responses, induced by both whole proteome of SARS-CoV-2 and specific viral proteins. A host-virus interactome was inferred, applying an explorative algorithm (Random Walk with Restart, RWR) triggered by 28 proteins of SARS-CoV-2. The analysis of PPI allowed to estimate the distribution of SARS-CoV-2 proteins in the host cell. Interactome built around one single viral protein allowed to define a different response, underlining as ORF8 and ORF3a modulated cardiovascular diseases and pro-inflammatory pathways, respectively. Finally, the network-based approach highlighted a possible direct action of ORF3a and NS7b to enhancing Bradykinin Storm. All results were reported in the paper Messina et al. 2021 (DOI: 10.1038/s41419-021-03881-8).
Many functions performed for this paper were based on funtions of R packages RandomWalkRestartMH (DOI: 10.1093/bioinformatics/bty637)  

Moreover, here two R scripts are reported with funtictions for interactome building and Gene Enrichment analysis, while a PPI database is reported.  
