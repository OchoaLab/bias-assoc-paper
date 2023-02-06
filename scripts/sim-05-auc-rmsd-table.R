# gathers small tables from each rep into a single big one (for easy commits and downstream analysis)

library(optparse) # for terminal options
library(readr)    # to read tables
library(dplyr)    # for bind_rows

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
herit <- opt$herit

if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# go where the data is
setwd( '../data/' )
setwd( dir_out )

# tibble to grow
data <- NULL

# include additional level if heritability is non-default
dir_herit <- '' # required this way for default herit
if ( herit != 0.8 )
    dir_herit <- paste0( 'h-', herit, '/' )

# load pre-calculated AUCs and SRMSDs
for ( rep in 1 : n_rep ) {
    file_rep <- paste0( 'rep-', rep, '/', dir_herit, 'eval.txt.gz' )
    data_rep <- read_tsv( file_rep, col_types = 'cid' )
    data <- bind_rows( data, data_rep )
}

# store final data in base, with herit sublevel if non-default
if ( dir_herit != '' ) {
    # create first time if needed
    if ( !dir.exists( dir_herit ) )
        dir.create( dir_herit )
    # move in there
    setwd( dir_herit )
}
write_tsv( data, 'eval.txt.gz' )
