#install.packages("rngtools", repos="http://R-Forge.R-project.org")
#install.packages("doRNG", repos="http://R-Forge.R-project.org")
#install.packages(""WebGestaltR")

# http://www.webgestalt.org/WebGestalt_2019_Manual.pdf

library("WebGestaltR")

testFunction <- function (...) {
  return(tryCatch(WebGestaltR(...), error=function(e) NULL))
}

#args = commandArgs(trailingOnly=TRUE)
#projectName =  args[1] # "cVamp_intersection-strict_fc" ### prefisso del file
#folderName =  args[2] #"../rh4_files/" ### directory principale che ha dentro webgestalt e GO
#program = args[3] ### nome directory di edger o deseq2


#print the list of organisms
#listOrganism(hostName = "http://www.webgestalt.org/")

#print the list of id types
#listIdType(organism = "hsapiens",hostName = "http://www.webgestalt.org/")

#print the list of gene sets
#listGeneSet(organism = "mmusculus",hostName = "http://www.webgestalt.org/")

ORA_function <- function(folderName,program,projectName,refList, geneList){

categories = c("geneontology_Biological_Process",
               "geneontology_Biological_Process_noRedundant",
"geneontology_Cellular_Component",
"geneontology_Cellular_Component_noRedundant",
"geneontology_Molecular_Function",
"geneontology_Molecular_Function_noRedundant",
#"pathway_KEGG",
#"pathway_Panther",
"pathway_Reactome",
#"network_CORUM",
"pathway_Wikipathway",
#"pathway_Wikipathway_cancer",
#"network_Kinase_phosphosite",
#"network_Kinase_target",
#"network_PPI_BIOGRID",
#"network_PTMsigDB",
"network_Transcription_Factor_target",
 "network_miRNA_target"#,
#"phenotype_Mammalian_Phenotype_Ontology",
#"phenotype_Human_Phenotype_Ontology",
#"chromosomalLocation_CytogeneticBand",
#"community-contributed_5htGeneSets_Conte",
#"community-contributed_MuscleGeneSets_Duddy_2017"
)
if (!dir.exists(paste0(folderName,program))) {
  dir.create(paste0(folderName,program))                  
}

dir.create(paste0(folderName,program,"/webgestalt/")) 
dir.create(paste0(folderName,program,"/webgestalt/ORA/"))

out.dir = mapply(paste, folderName,program,"/webgestalt/ORA/",categories,"/", MoreArgs=list(sep="") , SIMPLIFY=TRUE)
mapply(dir.create, out.dir)


### load the list of interesting genes and the background into a vector
#refList <- scan(paste(folderName,program,"/GO/",projectName,"_background.txt",sep=""), character(), quote = "")
#geneList <- scan(paste(folderName,program,"/GO/",projectName,".txt",sep=""), character(), quote = "")

wg = mapply(testFunction, enrichDatabase = categories , outputDirectory = out.dir, 
            MoreArgs=list(enrichMethod="ORA", organism="hsapiens",reportNum = 50,
                          interestGene=geneList, interestGeneType="genesymbol",
                          referenceGene=refList,referenceGeneType="genesymbol", 
                          minNum = 5, maxNum=2000, isOutput=TRUE, projectName=projectName) )

sink(paste(folderName,program,"/webgestalt/ORA/", projectName,".txt", sep="") )
print(wg)
sink()



# geneList <- scan(paste(folderName,program,"/GO/",projectName,"_up.txt",sep=""), character(), quote = "")
# 
# wg = mapply(testFunction, enrichDatabase = categories , outputDirectory = out.dir, MoreArgs=list(enrichMethod="ORA", organism="hsapiens",interestGene=geneList, interestGeneType="genesymbol",referenceGene=refList,referenceGeneType="ensembl_gene_id", minNum = 5, maxNum=2000, isOutput=TRUE, projectName=paste(projectName,"_up",sep="")) )
# 
# sink(paste(folderName,program,"/webgestalt/ORA/", projectName,"_up.txt", sep="") )
# print(wg)
# sink()
# 
# 
# geneList <- scan(paste(folderName,program,"/GO/",projectName,"_down.txt",sep=""), character(), quote = "")
# 
# wg = mapply(testFunction, enrichDatabase = categories , outputDirectory = out.dir, MoreArgs=list(enrichMethod="ORA", organism="hsapiens",interestGene=geneList, interestGeneType="genesymbol",referenceGene=refList,referenceGeneType="ensembl_gene_id", minNum = 5, maxNum=2000, isOutput=TRUE, projectName=paste(projectName,"_down",sep="")) )
# 
# sink(paste(folderName,program,"/webgestalt/ORA/", projectName,"_down.txt", sep="") )
# print(wg)
# sink()
}
