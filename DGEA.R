#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : combined_DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Differential Gene Expression Analysis
# ---------------------------------------------------------

#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

# source('http://bioconductor.org/biocLite.R')
# biocLite('GEOquery')
# install.packages('rjson')

#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

library('Cairo')
library('dendextend')
library('GEOquery')
library('ggplot2')
library('gplots')
library('limma')
library('plyr')
library('RColorBrewer')
library('reshape2')

# Command Line Arguments
working.dir   <- '/Users/nazrathnawaz/Dropbox/GroupPoject/RCore/CourseWorkScripts'
accession.id  <- 'GDS5093' #  GDS5092 GDS5091 GDS5088 GDS5086 GDS3795
factor_type   <- 'disease.state'
#factor_type  <- 'genotype/variation'
#factor_type  <- 'development.stage'
#factor_type  <- 'infection'


# import data sets and process into expression data
gse              <- getGEO(accession.id, GSEMatrix = TRUE)       # Load GEO data
eset             <- GDS2eSet(gse, do.log2=TRUE)                  # Convert into ExpressionSet Object
X                <- exprs(eset)                                  # Get Expression Data
geneNames        <- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
names(geneNames) <- rownames(X)

# Create a data frame with the factors
expression.info <- data.frame(pData(eset)[factor_type],
                              row.names     = sampleNames(eset))

data <- within(melt(X), {
    Factors = expression.info[Var2, factor_type]
})

# BoxPlot with colours based on the factors given
boxplot <- ggplot(data) + labs(x = 'Sample', y = 'Expression Levels') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'bottom'
    ) + geom_boxplot(aes(x = Var2, y = value, colour = Factors))

# store Boxplot as an .svg file in the working directory
filename <- paste(working.dir,"boxplot.svg",sep = "")
ggsave(filename, plot=boxplot, width = 8, height = 4)


######
# Differential Gene Expression
######
#which((pData(eset) == levels(pData(eset))[1]) == TRUE)

pClass        <- as.factor(pData(eset)[,factor_type])
names(pClass) <- sampleNames(eset)

# remove white spaces in factors and rename
levels(pClass) <-unlist(lapply(levels(pClass), function (x) gsub(" ", "", x, fixed = TRUE)))

arr   <- levels(pClass)
phase <- c()

## to be replaced ##
for(i in seq(1,length(arr)-1))
{
    for(j in seq(i+1,length(arr)))
    {
        a <- paste("pClass", arr[i],sep="")
        b <- paste("pClass", arr[j],sep="")
        result <- paste(a, b, sep="-")
        phase <- append(phase, result)
    }
}


#### Using the limma package ####

design  <- model.matrix(~0+pClass)

# plots linear model for each gene and estimate fold changes and standard errors
fit     <- lmFit(X, design)

# set contrasts for all classes
contrasts <- makeContrasts(contrasts=phase,
                           levels = design)

fit <- contrasts.fit(fit, contrasts)  

# empirical Bayes smoothing to standard errors
fit <- eBayes(fit)

## sort by p-values (to be reviewed) ##
toptable <- toptable(fit, sort.by='logFC', number=250, genelist = geneNames)
toptable['gene'] <- geneNames[rownames(toptable)]

# Create Sub Data
X.toptable <- X[rownames(toptable),]
rownames(X.toptable) <- geneNames[rownames(toptable)]

######
# Draw Graphs for DGEA
######

# store Heatmap as an .svg file in the working directory
filename <- paste(working.dir,"Heatmap.svg",sep = "")
CairoSVG(file = filename)
color_scale <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
heatmap1 <- heatmap.2(X.toptable, col=color_scale, scale='row', key=T, keysize=1.5,
                      dendrogram='row', density.info='none', trace='none', cexCol=0.5, cexRow=0.1)
dev.off()

# store Histogram as an .svg file in the working directory
filename <- paste(working.dir,"Histogram.svg",sep = "")
CairoSVG(file = filename)
hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', xlab = "Adjusted p-values")
dev.off()

# store Volcano plot as an .svg file in the working directory
filename <- paste(working.dir,"Volcano.svg",sep = "")
CairoSVG(file = filename)
volcanoplot(fit, coef=1, highlight=20, names=geneNames, col='steelblue', xlab='Log Fold Change',
            ylab='Log Odds', pch=16, cex=0.5)
dev.off()

# store Top20 genes as a .csv file in the working directory
filename <- paste(working.dir,"TopGenes.csv",sep = "")
write.csv(toptable[1:20,], file = filename)


