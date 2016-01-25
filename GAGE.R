#Analysis specific to the dengue dataset.



#---------------------------------

###Loading the data

#---------------------------------


source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(GEOquery)


#Importing data from GEO
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

#Get dataset with expression info
X <- Table(gse)

#Converting GSE to an expression set object
eset <- GDS2eSet(gse, do.log2=TRUE)


pDat <- pData(eset)

#------------------------------------

###Actually using GAGE

#------------------------------------


##Loading gage and associated gene sets
biocLite(c("gage","gageData","GO.db", "pathview" ))
library(gage) #Does the analysis
library(gageData) #Lets data be used by gage
library(pathview) #Visualises interaction networks & used to get ENTREZ IDs
library(GO.db) ##Downloads GO datasets


#---------------------
###Data prep
#---------------------


##Remove probe ID column & convert into data matrix
X1<- X
X1<- X[,-1]
X1_matrix<-data.matrix(X1)


##Create two column table containing entrez IDs for geodataset
id.map.refseq <- id2eg(ids = X$IDENTIFIER, category = "SYMBOL", org = "hsa")
#data(bods) - contains values  for 'org' argument. 


##Replace gene symbols with ENTREZ ID in dataset matrix
for (i in 1:length(id.map.refseq[,1])){
  if (id.map.refseq[i,1] == X1_matrix[i,1]){
   X1_matrix[i,1]<-id.map.refseq[i,2]
  }
}

##Remove rows without ENTREZ IDs
X1_matrix<-X1_matrix[complete.cases(X1_matrix),]

##Make first column rownames
GEOdataset <- X1_matrix[,-1]
rownames(GEOdataset) <- X1_matrix[,1]
##Convert to numerical matrix (for gage function)
class(GEOdataset) <- "numeric"  

#-------------------

##Carrying out GAGE 
#-------------------


#Get positions of specific, in vector form
cn=colnames(X1_matrix2)

##Find out which group each sample is part of
##Use this to form two groups of samples for further analysis


Group1<-c()
Group2<-c()

for (a in 1:length(pDat$sample)){
  if (pDat$infection[a] == "Dengue virus"){
    Group1_pos<- grep(pDat$sample[a], cn)
    Group1<-c(Group1, Group1_pos)
  }
  if (pDat$infection[a] == "control"){
    Group2<- c(Group2, (grep(pDat$sample[a], cn)))
  }
}

##Load gene sets
data(kegg.gs)
data(go.gs)


##Loading human kegg gene sets

kg.hsa=kegg.gsets("hsa") #this picks out the human sets

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step

save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object

##Loading human GO gene sets

go.hs=go.gsets(species="human")
go.bp=go.hs$go.sets[go.hs$go.subs$BP]
go.mf=go.hs$go.sets[go.hs$go.subs$MF]
go.cc=go.hs$go.sets[go.hs$go.subs$CC]
save(go.bp, go.mf, go.cc, file="go.hs.gsets.RData")
#for Bioconductor species supported by go.gsets function:
#data(bods)
#print(bods)



##Using the gage function to carry out analysis


#name           <- gage(data    ,  genesets used, control group, experimental group)
GEOdataset.kegg.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, compare= 'unpaired')

GEOdataset.go.p <- gage(GEOdataset, gsets = go.gs, ref = Group2, samp = Group1,compare= 'unpaired')

GEOdataset.bp.p <- gage(GEOdataset, gsets = go.bp,ref = Group2, samp = Group1, compare= 'unpaired')


#here we can look at the data
str(GEOdataset.kegg.p, strict.width='wrap')

head(GEOdataset.kegg.p$greater[, 1:5], 4)

head(GEOdataset.kegg.p$less[, 1:5], 4)

head(GEOdataset.kegg.p$stats[, 1:5], 4)


