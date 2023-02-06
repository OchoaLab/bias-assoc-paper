library(optparse) # for terminal options
library(readr)    # to read tables
library(tibble)
library(simtrait) # pval_aucpr

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "Replicate number", metavar = "int"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
rep <- opt$rep
herit <- opt$herit

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )
dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )
# enter herit subdir only if non-default
if ( herit != 0.8 ) {
    dir_herit <- paste0( 'h-', herit )
    setwd( dir_herit )
}

# load tibbles
pvals <- read_tsv( 'pvals.txt.gz', show_col_types = FALSE )

# and true causal info, for AUC
load( 'simtrait.RData' )

# calculate AUCs!
aucs <- lapply( pvals, function(x) pval_aucpr( x, causal_indexes ) )
# calculate RMSD_p
rmsds <- lapply( pvals, function(x) pval_srmsd( x, causal_indexes ) )

# write data into a tibble, for merging later!
data <- tibble(
    kinship = names( aucs ),
    rep = rep,
    rmsd = unlist( rmsds ),
    auc = unlist( aucs )
)
write_tsv( data, 'eval.txt.gz' )
