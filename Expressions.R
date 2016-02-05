#!/usr/bin/Rsript
# ---------------------------------------------------------------------------
# Filename      : DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Retrieve individual gene expressions and convert to JSON
# Run           : Rscript Expressions.R --dbrdata /Users/sureshhewapathirana/Desktop/topexpr.rData --outputdir /Users/sureshhewapathirana/Desktop/expression.json --rowid LOC100288410
# ---------------------------------------------------------------------------

library('argparser')    # Argument passing
library('jsonlite')     # Convert R object to JSON format

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--dbrdata", help="Full file path of rData file") 
parser <- add_argument(parser, "--outputdir", help="The outout directory where json file to besaved")    
parser <- add_argument(parser, "--rowid", help="Row Id of the X matrix")    

# Gene Expression specific parameters


# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                          Loading Saved Dataset                            #
#############################################################################


if (file.exists(argv$dbrdata)){
    load(file = argv$dbrdata)
}else{
    print("ERROR:File not found")
    q(save = "default")
}

#############################################################################
#                          Retrieve Expression data                         #
#############################################################################

if((!is.na(argv$outputdir))&&(!is.na(X.toptable))){
    index.group1 <- which((expression.info['population']== 'Group1')==TRUE)
    g1 <- list( x = names(X.toptable[argv$rowid, index.group1]),
                y = as.double(X.toptable[argv$rowid, index.group1]))
  
    index.group2 <- which((expression.info['population']== 'Group2')==TRUE)
    g2 <- list(x = names(X.toptable[argv$rowid, index.group2]),
                   y = as.double(X.toptable[argv$rowid, index.group2]))
    write(toJSON(list(list(group1 = g1 ,group2 = g2))), argv$outputdir)
}