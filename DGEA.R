#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : DGEA.R
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

#############################################################################
#                        Command Line Arguments                          #
#############################################################################

output.dir     <- "/Users/sureshhewapathirana/Desktop/"

# --------- Geo DataSet Input ------------

accession.id    <- 'GDS5093' #  GDS5092 GDS5091 GDS5088 GDS5086 GDS3795
factor.type     <- 'disease.state'
population1     <- c('Dengue Hemorrhagic Fever','Convalescent')
population2     <- c('healthy control')
pop.name1       <- "Dengue"
pop.name2       <- "Normal"
pop.colour1     <- "#b71c1c" # Red  
pop.colour2     <- "#0d47a1" # Blue 

# --------- Volcano Plot ------------
no.of.top.genes <- 250 
toptable.sortby <- "p" # or "logFC"

#############################################################################
#                        Testing Variables                         #
#############################################################################

#factor_type  <- 'genotype/variation'
#factor_type  <- 'development.stage'
#factor_type  <- 'infection'

#############################################################################
#                        GEO Input                          #
#############################################################################

# import data sets and process into expression data
gse              <- getGEO(accession.id, GSEMatrix = TRUE)       # Load GEO data
eset             <- GDS2eSet(gse, do.log2=TRUE)                  # Convert into ExpressionSet Object
X                <- exprs(eset)                                  # Get Expression Data
gene.names        <- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
names(gene.names) <- rownames(X)
pClass           <- pData(eset)[factor.type]
colnames(pClass) <- 'factor.type'
samples          <- rownames(pClass)



# Create a data frame with the factors
expression.info  <- data.frame(pClass, Sample = samples, row.names = samples)


expression.info <- within(expression.info, {
  Population        = ifelse (factor.type %in% population1,
                       'Group1', # if true
                       ifelse( factor.type %in% population2, 'Group2', NA) ) # if false
  population.colour = ifelse (factor.type %in% population1,
                       pop.colour1, # if true
                       ifelse( factor.type %in% population2, pop.colour2, '#000000') ) # if false
})
expression.info$Population <- as.factor(expression.info$Population)

data <- within(melt(X), {
    Phenotypes = expression.info[Var2, 'factor.type']
})

boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = Phenotypes)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = as.vector(expression.info$population.colour)), legend.position = 'right')+ labs(x = 'Samples', y = 'Expression Levels')
# store Boxplot as an .svg file in the working directory
filename <- paste(output.dir,"boxplot.svg",sep = "")
ggsave(filename, plot=boxplot, width = 8, height = 4)


######
# Differential Gene Expression
######


#pClass        <- as.factor(pData(eset)[,factor_type])
#names(pClass) <- sampleNames(eset)

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info),]
X <- X[,(colnames(X) %in% rownames(expression.info))]

# Created a new Phenotype class
newPClass           <- expression.info$Population
names(newPClass)    <- expression.info$Sample
newPClass

design  <- model.matrix(~0+newPClass)

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

toptable['gene'] <- gene.names[rownames(toptable)]

#Create Sub Data
X.toptable <- X[rownames(toptable),]
rownames(X.toptable) <- gene.names[rownames(toptable)]

######
# Draw Graphs for DGEA
######

# store Heatmap as an .svg file in the working directory
filename <- paste(output.dir,"Heatmap.svg",sep = "")
CairoSVG(file = filename)
color_scale <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
heatmap1 <- heatmap.2(X.toptable, col=color_scale, scale='row', key=T, keysize=1,
                      dendrogram='column', density.info='none', 
                      trace='none', cexCol=0.5, cexRow=0.1,
                      ColSideColors =  expression.info$population.colour,
                      Colv = TRUE, Rowv = TRUE)
dev.off()

# store Histogram as an .svg file in the working directory
filename <- paste(output.dir,"Histogram.svg",sep = "")
CairoSVG(file = filename)
hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', xlab = "Adjusted p-values", main=NULL)
dev.off()

# store Volcano plot as an .svg file in the working directory
filename <- paste(output.dir,"Volcano.svg",sep = "")
CairoSVG(file = filename)
with(toptable, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-1,1)))
#volcanoplot(fit, coef=1, highlight=20, names=gene.names, col='steelblue', xlab='Log Fold Change',
#            ylab='Log Odds', pch=16, cex=0.5)
dev.off()

# store Top20 genes as a .csv file in the working directory
filename <- paste(output.dir,"TopGenes.csv",sep = "")
write.csv(toptable[1:no.of.top.genes,], file = filename)

