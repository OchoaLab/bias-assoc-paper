library(optparse) # for terminal options
library(genio)    # to read GRM files
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF

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
setwd( 'kinship' )

# they have these names
# all are 2x, scaled as GCTA wants them, halve here for plots
kinship_true <- read_grm( 'true' )$kinship / 2
kinship_popkin <- read_grm( 'popkin' )$kinship / 2
kinship_std_rom <- read_grm( 'std_rom' )$kinship / 2
kinship_std_rom_lim <- read_grm( 'std_rom_lim' )$kinship / 2
kinship_std_mor <- read_grm( 'std_mor' )$kinship / 2
kinship_gcta <- read_grm( 'gcta' )$kinship / 2
kinship_gcta_lim <- read_grm( 'gcta_lim' )$kinship / 2
kinship_wg <- read_grm( 'wg' )$kinship / 2
kinship_wg_lim <- read_grm( 'wg_lim' )$kinship / 2

# save figure in lower level
setwd( '..' )

# visualize all matrices for test
dims <- fig_scale( ratio = 3/4 )
fig_start(
    'kinship',
    width = dims[1],
    height = dims[2]
)
plot_popkin(
    inbr_diag(
        list(
            kinship_true,
            kinship_popkin,
            NULL,
            kinship_std_rom_lim,
            kinship_std_rom,
            kinship_std_mor,
            kinship_gcta_lim,
            kinship_gcta,
            NULL,
            kinship_wg_lim,
            kinship_wg
        )
    ),
    titles = c(
        'Truth',
        'Popkin est.',
        'Standard ROM lim.',
        'Standard ROM est.',
        'Standard MOR est.',
        'GCTA lim.',
        'GCTA est.',
        'Weir-Goudet lim.',
        'Weir-Goudet est.'
    ),
    layout_rows = 4,
    mar = c(0, 2),
    panel_letters_adj = 0 # old default, works better here because there's no labels
)
fig_end()
