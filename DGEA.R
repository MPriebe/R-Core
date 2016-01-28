#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R 								   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa #
# Description   : Differential Gene Expression Analysis    # 
# Rscript DGEA.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005 --outputdir "/Users/sureshhewapathirana/Desktop/"
# ---------------------------------------------------------#



#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

# source('http://bioconductor.org/biocLite.R')
# biocLite('GEOquery')
# install.packages("argparser")
# install.packages('Cairo')
# install.packages('dendextend')
# install.packages('GEOquery')
# install.packages('ggplot2')
# install.packages('gplots')
# install.packages('rjson')
# install.packages('limma')
# install.packages('plyr')
# install.packages('RColorBrewer')
# install.packages('reshape2')

#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

library('argparser')
library('Cairo')
library('dendextend')
library('GEOquery')
library('ggplot2')
library('gplots')
library('rjson')
library('limma')
library('plyr')
library('RColorBrewer')
library('reshape2')

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")
parser <- add_argument(parser, "--accession", help="Accession Number of the GEO Database")
parser <- add_argument(parser, "--factor"     , help="input file")    # Factor type to be classified by
parser <- add_argument(parser, "--popA", nargs='+', help="input file")    # GroupA - all the selected phenotypes (atleast one)
parser <- add_argument(parser, "--popB", nargs='+', help="input file")    # GroupB - all the selected phenotypes (atleast one)
parser <- add_argument(parser, "--popname1"     , help="input file")    # name for GroupA
parser <- add_argument(parser, "--popname2"     , help="input file")    # name for GroupB
parser <- add_argument(parser, "--topgenecount"   , help="input file")    # number of top genes to be used
parser <- add_argument(parser, "--foldchange"   , help="input file")    # fold change cut off
parser <- add_argument(parser, "--thresholdvalue" , help="input file")    # threshold value cut off
parser <- add_argument(parser, "--outputdir"    , help="input file")    # GEO Accession ID
parser <- add_argument(parser, "--dbrdata"    , help="input file")    # GEO Accession ID

# allow arguments to be run via the command line
argv <- parse_args(parser)
output.dir     <- argv$outputdir

# --------- Geo DataSet Input ------------ #
accession.id    <- argv$accession           
factor.type     <- argv$factor           
population1     <- unlist(strsplit(argv$popA, ",")) 
population2     <- unlist(strsplit(argv$popB, ",")) 
pop.name1       <- argv$popname1         
pop.name2       <- argv$popname2        
pop.colour1     <- "#b71c1c"            # Red  
pop.colour2     <- "#0d47a1"            # Blue 


# --------- Volcano Plot ------------ #
no.of.top.genes <- as.numeric(argv$topgenecount)   # 250
toptable.sortby <- "p"                 				# sort by p-value (default)
fold.change 	<- as.numeric(argv$foldchange)       # 0.3
threshold.value <- as.numeric(argv$thresholdvalue)   # 0.005 # 0.05/no.of.top.genes -  Bonferroni cut-off



#if file.exists(argv$dbrdata){
load(file = argv$dbrdata)
#}else{
   # print("ERROR:File not found")
#}

#############################################################################
#                       Factor Selection                                 #
#############################################################################

gene.names        <- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
names(gene.names) <- rownames(X)
pClass            <- pData(eset)[factor.type]
colnames(pClass)  <- 'factor.type'
samples           <- rownames(pClass)

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Create a data frame with the factors
expression.info  <- data.frame(pClass, Sample = samples, row.names = samples)


# Introduce two columns to expression.info - 
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
expression.info <- within(expression.info, {
    population        = ifelse (factor.type %in% population1,'Group1', # if true
                                ifelse( factor.type %in% population2, 'Group2', NA) ) # if false
    population.colour = ifelse (factor.type %in% population1,pop.colour1, # if true
                                ifelse( factor.type %in% population2, pop.colour2, '#000000') ) # if false
})


# Convert to a factor
expression.info$population <- as.factor(expression.info$population)

data <- within(melt(X), {
    phenotypes = expression.info[Var2, 'factor.type']
})


# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info),]
X <- X[,(colnames(X) %in% rownames(expression.info))]

# Created a new Phenotype class
newPClass           <- expression.info$population
names(newPClass)    <- expression.info$Sample

#############################################################################
#                        Top Table                                    #
#############################################################################

find.toptable <- function(X, newPClass, toptable.sortby, no.of.top.genes, gene.names){
    
    design  <- model.matrix(~0 + newPClass)
    
    # plots linear model for each gene and estimate fold changes and standard errors
    fit     <- lmFit(X, design)
    
    # set contrasts for all classes
    contrasts <- makeContrasts(contrasts="newPClassGroup1-newPClassGroup2",
                               levels = design)
    
    fit <- contrasts.fit(fit, contrasts)
    
    # empirical Bayes smoothing to standard errors
    fit <- eBayes(fit)
    
    # Sort.by shoudl be variable - 'p' or 'LogFC'
    toptable <- topTable(fit, sort.by= toptable.sortby, number=no.of.top.genes, genelist = gene.names)
    
    return(toptable)
}


filtered.toptable <- function(toptable, gene.names){
    toptable['gene'] <- gene.names[rownames(toptable)]
    
    #Create Sub Data
    X.toptable <- X[rownames(toptable),]
    rownames(X.toptable) <- gene.names[rownames(toptable)]
    
    return(X.toptable)
}


#############################################################################
#                        Graphical Representations                          #
#############################################################################

# Initial Boxplot
samples.boxplot <- function(){
    boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = phenotypes)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = as.vector(expression.info$population.colour)), legend.position = 'right')+ labs(x = 'Samples', y = 'Expression Levels')
    # store Boxplot as an .png file in the working directory
    filename <- paste(output.dir,"boxplot.png",sep = "")
    ggsave(filename, plot=boxplot, width = 8, height = 4)
}

# Heatmap
heatmap <- function(X, sample.colours, cv = TRUE, rv = TRUE){
    # store Heatmap as an .png file in the working directory
    filename <- paste(output.dir,"heatmap.png",sep = "")
    CairoPNG(file = filename, width = 800, height = 800, pointsize = 12)
    color_scale <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
    heatmap1 <- heatmap.2(X, col=color_scale, scale='row', 
                          key=T, keysize=1,
                          dendrogram='column', density.info='none', 
                          trace='none', cexCol=0.6, cexRow=0.1,
                          ColSideColors = sample.colours,
                          Colv = cv, Rowv = rv)
    dev.off()
}

# Adjusted p-value barplot
adj.p.val.histogram <- function(toptable){
    # store Histogram as an .png file in the working directory
    filename <- paste(output.dir,"histogram.png",sep = "")
    CairoPNG(file = filename, width = 600, height = 600)
    hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', xlab = "Adjusted p-values", main=NULL)
    dev.off()
}

                                                #Bonferroni cut-off    
volcanoplot2 <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    
    # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    toptable$threshold = as.factor(abs(toptable$logFC) > fold.change & toptable$P.Value < t)
    
    # Construct the plot object
    vol = ggplot(data=toptable, aes(x=toptable$logFC, y=-log10(toptable$P.Value), colour=threshold)) +
        geom_point(alpha=0.4, size=1.75)  + xlim(c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1)) + ylim(c(0, max(-log10(toptable$P.Value))+0.5)) +
        xlab("log2 fold change") + ylab("-log10 p-value")
    filename <- paste(output.dir,"volcano.png",sep = "")
    ggsave(filename, plot=vol, height = 6, width = 6)
}

get.vol.data <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    
    vol <- data.frame(cbind(round(toptable$logFC,3), round(-log10(toptable$P.Value),3)))
    vol.list <- list(logFC = round(toptable$logFC,3),
                     pVal  = round(-log10(toptable$P.Value),3))
    return(vol.list)
}

# Top genes table
top.genes <- function(toptable, n){
    # store Top20 genes as a .csv file in the working directory
    filename <- paste(output.dir,"topgenes.csv",sep = "")
    write.csv(toptable[1:n,], file = filename)
}

# Clustering dendogram
clustering <- function(dist.method = "euclidean", clust.method = "average"){
    hc <- hclust(dist(t(X),dist.method), clust.method) 
    dend <- as.dendrogram(hc)
    labels_colors(dend) <- expression.info$population.colour[order.dendrogram(dend)]
    filename <- paste(output.dir,"cluster.png",sep = "")
    CairoPNG(file = filename, width = 800, height = 800, pointsize = 12)
    plot(dend, main = "Cluster Dendrogram", xlab = "Samples")
    dev.off()
}


# Principal Component Analysis

get.pc.data <- function(X){
    Xpca <- prcomp(t(X), scale= TRUE)
    s <- summary(Xpca)
    
    # Individual contribution of each princle component
    expVar <- s$importance[2,] * 100   # convert to %
    
    # Cumulative Variance
    cumVar <- s$importance[3,] * 100
    
    pcnames <-names(expVar)
    
    names(expVar) <- NULL
    names(cumVar) <- NULL
    
    results <- list(pcnames = pcnames,
                    expVar = expVar,
                    cumVar = cumVar)
    
    return(results)
}

#############################################################################
#                        Function Calling             		                #
#############################################################################

json.list <- list()
analysis.list <- c("Boxplot","Toptable","Volcano", "PCA","Heatmap", "Clustering")

if ("Boxplot" %in% analysis.list){
    samples.boxplot()
}

if ("Toptable" %in% analysis.list){
    toptable <- find.toptable(X, newPClass, toptable.sortby, no.of.top.genes, gene.names)
    X.toptable <- filtered.toptable(toptable, gene.names)
    adj.p.val.histogram(toptable)
    toptable.all <- find.toptable(X, newPClass, toptable.sortby, length(gene.names) , gene.names)
    json.list<- append(json.list,list(topgenes = toptable))
}

if ("Volcano" %in% analysis.list){
    volcanoplot2(toptable.all,fold.change)
    volcanoplot.data <- get.vol.data(toptable,fold.change)
    json.list<- append(json.list, list(vol = volcanoplot.data))
}

if ("PCA" %in% analysis.list){
    pcdata <- get.pc.data(X)
    json.list<- append(json.list, list(pc = pcdata))
}

if ("Heatmap" %in% analysis.list){
    heatmap(X.toptable,expression.info$population.colour)
}

if ("Clustering" %in% analysis.list){
    clustering()
}

if(length(json.list) != 0){
    filename <- paste(output.dir,"data.json",sep = "")
    write(toJSON(json.list), filename)
}
