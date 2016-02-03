#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : Installations.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Pre-Installation libraries
# ---------------------------------------------------------

#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

source('http://bioconductor.org/biocLite.R')
biocLite('GEOquery')
install.packages('argparser')
install.packages('Cairo')
install.packages('dendextend')
install.packages('GEOquery')
install.packages('ggplot2')
install.packages('gplots')
install.packages('heatmap3')
install.packages('rjson')
install.packages('limma')
install.packages('plyr')
install.packages('RColorBrewer')
install.packages('reshape2')