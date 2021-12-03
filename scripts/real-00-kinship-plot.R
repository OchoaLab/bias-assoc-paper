library(optparse) # for terminal options
library(genio)    # to read GRM files
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF
library(readr)

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )

# real data has the issue or sample ordering!
# read FAM for current order and subpop correspondence
fam <- read_fam( 'data' )
# read annotations
subpop_info <- read_tsv( 'pops-annot.txt', comment = '#', show_col_types = FALSE )
# map subpopulations using sub-subpopulations
fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]
# reorder individuals so subpopulations come out in desired order:
indexes <- order( match( fam$fam, subpop_info$pop ) )

setwd( 'kinship' )

# load kinship matrices
# all are 2x, scaled as GCTA wants them, halve here for plots
kinship_popkin <- read_grm( 'popkin' )$kinship / 2
kinship_std_rom <- read_grm( 'std_rom' )$kinship / 2
kinship_std_mor <- read_grm( 'std_mor' )$kinship / 2
kinship_gcta <- read_grm( 'gcta' )$kinship / 2
kinship_wg <- read_grm( 'wg' )$kinship / 2

# apply reordering to everything!
fam <- fam[ indexes, ]
kinship_popkin <- kinship_popkin[ indexes, indexes ]
kinship_std_rom <- kinship_std_rom[ indexes, indexes ]
kinship_std_mor <- kinship_std_mor[ indexes, indexes ]
kinship_gcta <- kinship_gcta[ indexes, indexes ]
kinship_wg <- kinship_wg[ indexes, indexes ]

# gather all the data, apply diagonal transformation
data <- inbr_diag(
    list(
        kinship_popkin,
        kinship_std_rom,
        kinship_std_mor,
        kinship_wg,
        kinship_gcta
    )
)

# the real kinship matrix generally has higher values, we want to preserve that range
# the biased kinship matrice have usually lower values but the diagonal has much larger values, those are the ones we want to cap
alpha <- 0.01
kinship_max <- quantile( diag( data[[1]] ), probs = 1 - alpha )
#kinship_max <- max( data[[1]] )
cap_kinship <- function( x ) {
    x[ x > kinship_max ] <- kinship_max
    return( x )
}
data <- lapply( data, cap_kinship )

# save figure in lower level
setwd( '..' )

# visualize all matrices for test
dims <- fig_scale( ratio = 3/4 )
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2] / 2 # hack
)
plot_popkin(
    data,
    titles = c(
        'Popkin',
        'Standard ROM',
        'Standard MOR',
        'Weir-Goudet',
        'GCTA'
    ),
    labs = fam$superpop,
    labs_line = 0.2,
    labs_cex = 0.7,
    layout_rows = 2,
    mar = c(2, 2)
#    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
