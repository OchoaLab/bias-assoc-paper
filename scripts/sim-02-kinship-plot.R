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
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option(c("-r", "--rep"), type = "integer", default = 1, 
                help = "Replicate number", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
rep <- opt$rep

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )

# go where the data is
setwd( '../data/' )
setwd( dir_out )
dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )
setwd( 'kinship' )

# load pre-existing data

# kinship files have these names, will appear in this order
codes <- c(
    'true',
    'popkin_rom',
    'popkin_mor',
    'std_rom_lim',
    'std_rom',
    'std_mor',
    'wg_rom_lim',
    'wg_rom',
    'wg_mor'
)
# map to human-readable names from table
titles <- kinship_methods$nice[ match( codes, kinship_methods$code ) ]

# read in all kinship matrices
data <- lapply( codes, function ( name ) {
    # all are 2x, scaled as GCTA wants them, halve here for plots
    read_grm( name )$kinship / 2
})

# save figure in lower level
setwd( '../..' )

# visualize all matrices for test
dims <- fig_scale( ratio = 1 )
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2]
)
plot_popkin(
    inbr_diag( data ),
    titles = titles,
    layout_rows = 3,
    mar = c(0, 2),
    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
