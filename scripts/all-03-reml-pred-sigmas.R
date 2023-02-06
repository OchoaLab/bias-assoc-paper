library(readr)
library(ochoalabtools)
library(popkin)
library(optparse)    # for terminal options

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
herit <- opt$herit

# include additional level if heritability is non-default
dir_herit <- '' # so stuff works for default case too
if ( herit != 0.8 )
    dir_herit <- paste0( 'h-', herit, '/' )
file_in <- paste0( dir_herit, 'reml.txt.gz' )

#################
### LOAD DATA ###
#################

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# go where the data is
# load precomputed sigmas
setwd( '../data/' )

# simulation first
# here's a weird hack to load old data for default (high) herit (so prev results are unchanged), but newer data for non-default herit (couldn't use old data for these new results)
if ( dir_herit != '' ) {
    setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-g20' )
} else {
    setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes' )
}
data_sim <- read_tsv( file_in, col_types = 'icdddddddd' )
setwd( '..' )

# now real data
setwd( 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01' )
data_real <- read_tsv( file_in, col_types = 'icdddddddd' )
setwd( '..' )

# store outputs in herit output
# create if it didn't already exist
if ( dir_herit != '' ) {
    if ( !dir.exists( dir_herit ) )
        dir.create( dir_herit )
    setwd( dir_herit )
}


# visualize raw sigma distributions
# these are useful because the simulations are replicates, so all the same values are being estimated over and over, and the variance is just noise in theory
boxplot_sigmas <- function ( data, main = '', leg = FALSE ) {
    # types vary per dataset, so determine them here
    types <- unique( data$type )
    # and by component/bias types (constants)
    components <- c( 'e', 'g' )
    short_names <- c( 'tr', 'std', 'wg' )
    sigmas <- list()
    labs_type <- c() # grow it in loop
    labs_short <- c() # ditto
    labs_component <- c() # ditto
    # this loop determines boxplot order
    for ( component in components ) {
        for ( type in types ) {
            for ( short_name in short_names ) {
                # name as it appears in column of `data`
                col_name <- paste0( component, '_', short_name )
                # separate data as desired
                sigmas[[ paste0( col_name, '_', type ) ]] <- data[[ col_name ]][ data$type == type ]
                # add pretty names for plot for these properties
                labs_component <- c( labs_component, if ( component == 'e' ) 'Residual' else 'Genetic' )
                labs_type <- c( labs_type, kinship_methods$type[ kinship_methods$code == paste0( 'std_', type ) ] )
                # this is how it matches my other tables
                if ( short_name == 'tr' ) {
                    short_name <- if ( type == 'rom_lim' ) 'true' else 'popkin_rom'
                } else
                    short_name <- paste0( short_name, '_rom' )
                labs_short <- c( labs_short, kinship_methods$short[ kinship_methods$code == short_name ] )
            }
        }
    }
    # make boxplot!
    boxplot(
        sigmas,
        main = main,
        names = NA, # individual labels will be plotted with rest
        xaxt = 'n',
        ylab = 'Variance estimate'
    )
    # reuse popkin-style labeling!  Great for hierarchical/factorial setup
    # though popkin has defaults for these, the raw function doesn't have defaults!
    popkin:::print_labels_multi(
                 labs = cbind( labs_short, labs_type, labs_component ),
                 labs_cex = c(0.7, 0.6, 1),
                 labs_las = c(2, 0, 0),
                 labs_line = c(0.5, 4, 6),
                 labs_lwd = 1, # default
                 labs_sep = c(FALSE, TRUE, TRUE),
                 labs_even = c(FALSE, TRUE, TRUE),
                 labs_ticks = FALSE, # default
                 labs_text = TRUE, # default
                 labs_col = 'black', # default
                 # these align with barplot-specific changes
                 xc_ind = 1 : length( sigmas ),
                 doMat = FALSE
             )
    # add some guide lines reflecting biases
    # this works because the true overal sigma is 1
    sigma_sq_e_true <- 1 - herit
    abline( h = herit, lty = 2, col = 'blue' )
    abline( h = sigma_sq_e_true, lty = 2, col = 'red' )
    if ( leg )
        legend(
            if ( herit > 0.5 ) 'topleft' else 'bottomleft',
            c('Genetic', 'Residual'),
            lty = 2,
            col = c('blue', 'red'),
            cex = 0.7,
            title = 'True variances'
        )
}

boxplot_errors <- function ( data, main = '' ) {
    # calculate errors
    # for environmental, prediction is that they are all the same
    data$err_e_std <- data$e_std - data$e_tr
    data$err_e_wg <- data$e_wg - data$e_tr
    data$err_e_both <- data$e_std - data$e_wg
    # for genetic variance, use factor predicted from theory
    data$err_g_std = data$g_std - data$g_tr * data$c_std
    data$err_g_wg = data$g_wg - data$g_tr * data$c_wg
    data$err_g_both = data$g_std - data$g_wg * ( data$c_std / data$c_wg )

    # in this plot only shorten GW and True
    kinship_methods$short[ kinship_methods$short== 'Weir-Goudet' ] <- 'WG'
    kinship_methods$short[ kinship_methods$short== 'True Kinship' ] <- 'True'
    
    # types vary per dataset, so determine them here
    types <- unique( data$type )
    # and by component/bias types (constants)
    components <- c( 'e', 'g' )
    short_names <- c( 'std', 'wg', 'both' )
    errors <- list()
    labs_type <- c() # grow it in loop
    labs_short <- c() # ditto
    labs_component <- c() # ditto
    # this loop determines boxplot order
    for ( component in components ) {
        for ( type in types ) {
            for ( short_name in short_names ) {
                # name as it appears in column of `data`
                col_name <- paste0( 'err_', component, '_', short_name )
                # separate data as desired
                errors[[ paste0( col_name, '_', type ) ]] <- data[[ col_name ]][ data$type == type ]
                # add pretty names for plot for these properties
                labs_component <- c( labs_component, if ( component == 'e' ) 'Residual' else 'Genetic' )
                labs_type <- c( labs_type, kinship_methods$type[ kinship_methods$code == paste0( 'std_', type ) ] )
                if ( short_name == 'both' ) {
                    lab_short <- 'Standard - WG'
                } else {
                    short_name_ref <- if ( type == 'rom_lim' ) 'true' else paste0( 'popkin_', type )
                    lab_short <- paste0(
                        kinship_methods$short[ kinship_methods$code == paste0( short_name, '_rom' ) ],
                        ' - ',
                        kinship_methods$short[ kinship_methods$code == short_name_ref ]
                    )
                }
                labs_short <- c( labs_short, lab_short )
            }
        }
    }
    # make boxplot!
    boxplot(
        errors,
        main = main,
        names = NA, # individual labels will be plotted with rest
        xaxt = 'n',
        ylab = 'Prediction error'
    )
    # reuse popkin-style labeling!  Great for hierarchical/factorial setup
    # though popkin has defaults for these, the raw function doesn't have defaults!
    popkin:::print_labels_multi(
                 labs = cbind( labs_short, labs_type, labs_component ),
                 labs_cex = c(0.7, 0.6, 1),
                 labs_las = c(2, 0, 0),
                 labs_line = c(0.5, 5.5, 7),
                 labs_lwd = 1, # default
                 labs_sep = c(FALSE, TRUE, TRUE),
                 labs_even = c(FALSE, TRUE, TRUE),
                 labs_ticks = FALSE, # default
                 labs_text = TRUE, # default
                 labs_col = 'black', # default
                 # these align with barplot-specific changes
                 xc_ind = 1 : length( errors ),
                 doMat = FALSE
             )
}

width <- fig_width()
fig_start(
    'preds-reml-sigmas',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 7
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_sigmas( data_sim, main = 'Admixed Family sim.', leg = TRUE )
panel_letter('A', adj = -0.2)
boxplot_sigmas( data_real, main = '1000 Genomes' )
panel_letter('B', adj = -0.2)
mtext( 'Variance Component, Kinship Estimate', side = 1, outer = TRUE )
fig_end()

# interpretation:
# - sigmas have a lot of noise, and some apparent biases by type, though the point of the paper is not really to characterize these biases, just to characterize their relative values across kinship bias types.  Also, boxplots show medians, not means, so bias not actually completely evident here.

fig_start(
    'preds-reml-errors',
    width = width,
    height = width / 2,
    mar_t = 1.5,
    mar_b = 8
)
par( mfrow = c(1,2) )
par( oma = c(1, 0, 0, 0) )
boxplot_errors( data_sim, main = 'Admixed Family sim.' )
panel_letter('A', adj = -0.2)
boxplot_errors( data_real, main = '1000 Genomes' )
panel_letter('B', adj = -0.2)
mtext( 'Variance Component, Kinship Estimate Pair', side = 1, outer = TRUE )
fig_end()

# interpretation
# - outliers are only real g ROM, and for both bias types the outliers are close, so popkin ROM must be the problem (proven in update including "Both" case)
#   - those outliers are very negative, meaning that popkin ROM's sigma is overestimated (too large) vs theory
