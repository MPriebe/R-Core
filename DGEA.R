#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript DGEA.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005 --xaxis "PC1" --yaxis "PC2" --distance "euclidean" --clustering "average" --dbrdata "~/Desktop/GDS5093.rData" --outputdir "/Users/sureshhewapathirana/Desktop/" --heatmaprows 100 --dendrow TRUE --dendcol TRUE --analyse "Boxplot,Volcano,PCA,Heatmap"
# ---------------------------------------------------------#

#############################################################################
#                        Import Libraries                                   #
#############################################################################

# silent library loading messages on command line
suppressMessages(library("limma"))
suppressMessages(library("dendextend"))
suppressMessages(library("gplots"))
suppressMessages(library("GEOquery"))

# load required libraries
library('argparser')    # Argument passing
library('Cairo')        # Plots saving
library('dendextend')   # Create dendogram
library('DMwR')         # Outlier Prediction for clustering
library('GEOquery')     # GEO dataset Retrieval
library('ggplot2')      # Graphs designing
library('gplots')       # Graphs designing
library('jsonlite')     # Convert R object to JSON format
library('pheatmap')     # Heatmap Generating
library('limma')        # Differencial Gene Expression Analysis
library('plyr')         # Splitting, Applying and Combining Data
library('RColorBrewer') # Import Colour Pallete
library('reshape2')     # Prepare dataset for ggplot


#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--outputdir", help="The outout directory where graphs get saved")
parser <- add_argument(parser, "--dbrdata", help="Downloaded GEO dataset full path")
parser <- add_argument(parser, "--analyse", help="List of analysis to be performed",nargs='+')

# Sample Parameters
parser <- add_argument(parser, "--accession", help="Accession Number of the GEO Database")
parser <- add_argument(parser, "--factor", help="Factor type to be classified by")   
parser <- add_argument(parser, "--popA", help="GroupA - all the selected phenotypes (atleast one)", nargs='+')
parser <- add_argument(parser, "--popB", help="GroupB - all the selected phenotypes (atleast one)", nargs='+')
parser <- add_argument(parser, "--popname1", help="name for GroupA")
parser <- add_argument(parser, "--popname2", help="name for GroupB")

# Volcano plot Parameters
parser <- add_argument(parser, "--foldchange", help="fold change cut off")
parser <- add_argument(parser, "--topgenecount", help="number of top genes to be used")
parser <- add_argument(parser, "--thresholdvalue" , help="threshold value cut off")

# Principle Component Parameters
parser <- add_argument(parser, "--xaxis", help="PC used as x axis")
parser <- add_argument(parser, "--yaxis", help="PC used as y axis")

# Clustering
parser <- add_argument(parser, "--distance", help="PC used as x axis")
parser <- add_argument(parser, "--clustering", help="PC used as y axis")

# Heatmap
parser <- add_argument(parser, "--heatmaprows", help="Number of genes show in the heatmap")
parser <- add_argument(parser, "--dendrow", help="Boolean value for display dendogram for Genes")
parser <- add_argument(parser, "--dendcol", help="Boolean value for display dendogram for Samples")

# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                        Command Line Arguments Retrieval                   #
#############################################################################

# General Parameters
output.dir      <- argv$outputdir
dbrdata         <- argv$dbrdata
analysis.list   <- unlist(strsplit(argv$analyse, ","))

# Sample Parameters
factor.type     <- argv$factor
population1     <- unlist(strsplit(argv$popA, ","))
population2     <- unlist(strsplit(argv$popB, ","))
pop.name1       <- argv$popname1
pop.name2       <- argv$popname2
pop.colour1     <- "#b71c1c"  # Red
pop.colour2     <- "#0d47a1"  # Blue

# Volcano plot Parameters
no.of.top.genes <- as.numeric(argv$topgenecount)
fold.change     <- as.numeric(argv$foldchange)
threshold.value <- as.numeric(argv$thresholdvalue)
toptable.sortby <- "p"

# Principle Component Parameters
x.axis <- argv$xaxis
y.axis <- argv$yaxis

# Clustering
if(argv$distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")){
    dist.method <- argv$distance
}else{
    dist.method <- "euclidean"
}

if(argv$clustering %in% c("ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid")){
    clust.method <- argv$clustering
}else{
    clust.method <- "average"
}

# Heatmap
heatmap.rows <- as.numeric(argv$heatmaprows)
cv <- as.logical(argv$dendrow)
rv <- as.logical(argv$dendcol)

# Remove command line argument variables
remove(parser)
remove(argv)

#############################################################################
#                        Testing                   #
#############################################################################

# factor.type     <- "disease.state"
# population1     <- c("Dengue Hemorrhagic Fever","Dengue Fever","Convalescent")
# population2     <- c("healthy control")
# 
# pop.name1       <- "Dengue"
# pop.name2       <- "Normal"
# pop.colour1     <- "#b71c1c"  # Red
# pop.colour2     <- "#0d47a1"  # Blue
# 
# # --------- Volcano Plot ------------ #
# no.of.top.genes <- 250 #as.numeric(argv$topgenecount)
# toptable.sortby <- "p"
# fold.change 	<- 0.3 #as.numeric(argv$foldchange)
# threshold.value <- 0.05 #as.numeric(argv$thresholdvalue)
# dbrdata         <-"/Users/sureshhewapathirana/Desktop/GDS5093.rData"
# output.dir      <-"/Users/sureshhewapathirana/Desktop/"
# x.axis <- "PC1"
# y.axis <- "PC2"
# dist.method <- "euclidean"
# clust.method <- "average"
# cv <- TRUE
# rv <- TRUE
# heatmap.rows <- 100
# analysis.list <- c("Boxplot","Volcano", "PCA","Clustering", "Heatmap")
#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################


if (file.exists(dbrdata)){
    load(file = dbrdata)
} else {
    if (is.null(argv$geodbpath)) {
        gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE) # Load data from downloaded file
    } else {
        gse <- getGEO(argv$accession, GSEMatrix = TRUE)            # Automatically Load GEO dataset
    }
    eset <- GDS2eSet(gse, do.log2=TRUE)                            # Convert into ExpressionSet Object
}

X    <- exprs(eset)                                            	 # Get Expression Data

## auto-detect if data needs transformation and log2 transform if needed
qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { X[which(X <= 0)] <- NaN
exprs(eset) <- log2(X) }


#############################################################################
#                       Factor Selection                                    #
#############################################################################

# Store gene names
X               <- exprs(eset)
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names

# Phenotype Selection
pClass          <- pData(eset)[factor.type]
colnames(pClass)<- 'factor.type'

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Create a data frame with the factors
expression.info  <- data.frame(pClass,
                               Sample = rownames(pClass),
                               row.names = rownames(pClass))

# Introduce two columns to expression.info -
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
expression.info <- within(expression.info, {
    population        = ifelse (factor.type %in% population1,'Group1', # if true
                                ifelse( factor.type %in% population2, 'Group2', NA) ) # if false
    population.colour = ifelse (factor.type %in% population1,pop.colour1, # if true
                                ifelse( factor.type %in% population2, pop.colour2, '#000000') ) # if false
})

# Convert population column to a factor
expression.info$population <- as.factor(expression.info$population)

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info),]
X <- X[,(colnames(X) %in% rownames(expression.info))]

# Data preparation for ggplot-Boxplot
data <- within(melt(X), {
    phenotypes = expression.info[Var2, 'factor.type']
    Groups = expression.info[Var2, 'population.colour']
})

# Created a new Phenotype class
newPClass           <- expression.info$population
names(newPClass)    <- expression.info$Sample

#############################################################################
#                        Top Table                                          #
#############################################################################

find.toptable <- function(X, newPClass, toptable.sortby, no.of.top.genes){
    
    design  <- model.matrix(~0 + newPClass)
    
    # plots linear model for each gene and estimate fold changes and standard errors
    fit     <- lmFit(X, design)
    
    # set contrasts for two groups
    contrasts <- makeContrasts(contrasts="newPClassGroup1-newPClassGroup2",
                               levels = design)
    
    fit <- contrasts.fit(fit, contrasts)
    
    # empirical Bayes smoothing to standard errors
    fit <- eBayes(fit)
    
    # Create top Table
    toptable <- topTable(fit,
                         sort.by= toptable.sortby, # 'p' or 'LogFC'
                         number=no.of.top.genes)
    
    return(toptable)
}

#############################################################################
#                        Graphical Representations                          #
#############################################################################

# Boxplot
samples.boxplot <- function(data){
    boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = Groups)) + theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = 'right')+ labs(x = 'Samples', y = 'Expression Levels') + scale_color_manual(name="Groups",values = c(pop.colour1,pop.colour2), labels=c(pop.name1,pop.name2))
    filename <- paste(output.dir,"boxplot.png",sep = "")
    ggsave(filename, plot=boxplot, width = 8, height = 4)
}

# Heatmap
heatmap <- function(X.matix, heatmap.rows = 100, cv = TRUE, rv = TRUE, dist.method, clust.method){
    X.matix<-X.toptable
    col.pal <- colorRampPalette(rev(brewer.pal(11, 'RdYlGn')))(100)
#     palette <- colorRampPalette(c("white","red"))(100)
#     outlier.colours <- palette[cut(outliers, 100)]
    
    # Column dendogram
    if(cv == TRUE){
        hc <- hclust(dist(t(X.matix), method = dist.method), method = clust.method)
        plot(as.dendrogram(hc))
        # Rank outliers using distance and clustering parameters
        o <- outliers.ranking(t(X.matix), 
                              test.data = NULL, 
                              method = "sizeDiff", # Outlier finding method
                              method.pars = NULL,
                              clus = list(dist = dist.method,
                                          alg = "hclust",
                                          meth = clust.method))
        
        annotation_col <- data.frame( Population = expression.info[,'population'],
                                      Factor = expression.info[,'factor.type'], 
                                      Outliers = o$prob.outliers)
        column.gap = 0
    }else{
        hc <- FALSE
        
        annotation_col <- data.frame( Population = expression.info[,'population'],
                                      Factor = expression.info[,'factor.type'])
        column.gap = length((which(annotation_col[,'Population'] == 'Group1')== TRUE))
       
    }
    
    rownames(annotation_col) = expression.info[,'Sample']
    #anno_colors = list(Outliers = outlier.colours)
    
    filename <- paste(output.dir,"heatmap.svg",sep = "")
    CairoSVG(file = filename)
    
    pheatmap(X.matix[1:heatmap.rows,], 
                       cluster_row = cv,
                       cluster_cols = hc,
                       annotation_col = annotation_col,
                       color = col.pal, 
                       fontsize = 6.5,
                       fontsize_row=3.5, 
                       fontsize_col = 3.5,
                       gaps_col=column.gap,
                       cutree_rows = 2)
    dev.off()
}

# Apply Bonferroni cut-off as the default thresold value
volcanoplot <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    
    # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    toptable$threshold = as.factor(abs(toptable$logFC) > fold.change & toptable$P.Value < t)
    
    # Construct the plot object
    vol = ggplot(data=toptable, aes(x=toptable$logFC, y=-log10(toptable$P.Value), colour=threshold)) +
        geom_point(alpha=0.4, size=1.75)  + xlim(c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1)) + ylim(c(0, max(-log10(toptable$P.Value))+0.5)) +
        xlab("log2 fold change") + ylab("-log10 p-value")
    
    # File saving
    filename <- paste(output.dir,"volcano.png",sep = "")
    ggsave(filename, plot=vol, height = 6, width = 6)
}

get.vol.data <- function(toptable){
    vol.list <- list( genes = toptable$ID,
                      logFC = round(toptable$logFC, 3),
                      pVal  = round(-log10(toptable$P.Value), 3) )
    return(vol.list)
}

# Principal Component Analysis
get.pc.data <- function(Xpca){
    
    s <- summary(Xpca)
    
    # Individual contribution of each princle component in percentages
    expVar <- s$importance[2,] * 100
    
    # Cumulative Variance in percentages
    cumVar <- s$importance[3,] * 100
    
    # PC names
    pcnames <-names(expVar)
    
    names(expVar) <- NULL
    names(cumVar) <- NULL
    
    results <- list(pcnames = pcnames,
                    expVar = expVar,
                    cumVar = cumVar)
    
    return(results)
}

get.pc.scatterplot.data <- function(Xpca){
    
    Xscores <- Xpca$x
    
    # Convert each column into list
   return(split(Xscores, rep(1:ncol(Xscores), each = nrow(Xscores))))
}

#############################################################################
#                        Function Calling                                 #
#############################################################################

json.list <- list()

# Toptable
toptable <- find.toptable(X, newPClass, toptable.sortby, no.of.top.genes)

# Adding to JSON file
temp.toptable <- toptable
names(temp.toptable) <- NULL
json.list<- append(json.list,list(tops = temp.toptable))

# Filter toptable data from X
X.toptable <-X[as.numeric(rownames(toptable)),]

if ("Boxplot" %in% analysis.list){
    samples.boxplot(data)
}

if ("Volcano" %in% analysis.list){
    
    # Create complete volcano plot
    toptable.all <- find.toptable(X, newPClass, toptable.sortby, length(gene.names))
    volcanoplot(toptable.all,fold.change)
    
    # save volcanoplot top data as JSON
    volcanoplot.data <- get.vol.data(toptable)
    json.list<- append(json.list, list(vol = volcanoplot.data))
}

if ("PCA" %in% analysis.list){
    
    Xpca <- prcomp(t(X.toptable), scale= TRUE)
    
    # PC individual and cumulative values
    pcdata <- get.pc.data(Xpca)
    json.list<- append(json.list, list(pc = pcdata))
    
    # PC scatter plot
    pcplotdata <- get.pc.scatterplot.data(Xpca)
    
    # adding both data to the json list
    json.list<- append(json.list, list(pcdata = pcplotdata))
}

if ("Heatmap" %in% analysis.list){
    heatmap(X.toptable, heatmap.rows = heatmap.rows, rv = rv, cv = TRUE, dist.method, clust.method)
}

if(length(json.list) != 0){
    output.dir <- "/Users/sureshhewapathirana/Desktop/"
    filename <- paste(output.dir,"data.json",sep = "")
    write(toJSON(json.list), filename)
}


