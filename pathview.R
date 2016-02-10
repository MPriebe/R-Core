#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : pathview.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : create interaction network    #
# Rscript Pathview.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.rData --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --rundir "~/Desktop/" --dev TRUE
# ---------------------------------------------------------#

# Load dependencies

library("argparser")    # Argument passing
library("pathview")     # Creates Heatmap



# Add parsers

parser <- add_argument(parser, "--gagedata", 
                       help="The GAGE results full path")

parser <- add_argument(parser, "--B.is.control",
                       help="Boolean for whether Group B should be treated as a control")

parser <- add_argument(parser, "--geneset", 
                       help="The KEGG geneset code")

parser <- add_argument(parser, "--rundir", 
                       help="The output directory where graphs get saved")




parser <- add_argument(parser, "--popA",
                       help = "GroupA - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--popB",
                       help = "GroupB - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--factor", 
                       help="Factor type to be classified by")



# allows arguments to be run via the command line
argv <- parse_args(parser)





