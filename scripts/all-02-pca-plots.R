library(genio)       # to write BED files for external software
library(ochoalabtools) # for nice PDF
library(readr)

# read in all kinship matrices
load_pcs <- function ( codes, K = 2 ) {
    data <- lapply( codes, function ( code ) {
        # all are 2x, scaled as GCTA wants them, halve here for plots
        kinship <- read_grm( code )$kinship / 2
        # calculate PCs, store only those
        pcs <- eigen( kinship )$vectors
        # subset further, we generally only want first two PCs
        # however, for sim and true/popkin the first PC is weird, so keep three for all instead
        pcs <- pcs[ , 1 : K ]
    })
    names( data ) <- codes
    return( data )
}

# make sure data is in the same order (as eigenvectors have arbitrary signs)
# apply to all together, align to the first one
align_pcs <- function( data, K = 2 ) {
    pcs1 <- data[[ 1 ]]
    for ( i in 2 : length( data ) ) {
        pcsi <- data[[ i ]]
        for ( j in 1 : K ) {
            # switch sign of PC if it is negatively correlated with that of first dataset
            if ( cor( pcs1[ , j ], pcsi[ , j ] ) < 0 )
                pcsi[ , j ] <- - pcsi[ , j ]
        }
        data[[ i ]] <- pcsi
    }
    return( data )
}

# plot PCs, several methods in a single panel
# `codes` helps plot subsets only (makes sense as we have triplets)
plot_pcs <- function( data, codes, xlim, ylim, colors = 1 : length( codes ), pch = '.', main = '', leg = FALSE ) {
    # set up area
    plot(
        NA,
        xlab = 'PC1',
        ylab = 'PC2',
        xlim = xlim,
        ylim = ylim
    )
    # for consistency, set main this way (shifts up a little bit)
    mtext( main, side = 3, line = 0.5 )
    # add points
    for ( code in codes ) {
        pcs <- data[[ code ]]
        i <- which( code == codes )
        points( pcs[ , 1 ], pcs[ , 2 ], col = colors[ i ], pch = pch )
    }
    # add legend
    if (leg ) 
        legend(
            'bottomright',
            kinship_methods$short[ match( codes, kinship_methods$code ) ],
            pch = NA,
            text.col = colors,
            cex = 0.8
        )
}

#################
### LOAD DATA ###
#################

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cc' )
# in simulation, all ROM limits were very similar to ROM estimates, so we'll just show estimates (only ones available in real data anyway)
kinship_methods <- kinship_methods[ kinship_methods$type != 'ROM lim.', ]
codes <- kinship_methods$code

# go where the data is
setwd( '../data/' )

# simulation first, use first replicate only
setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes/rep-1/kinship' )
message( 'Loading simulated data...' )
data_sim <- load_pcs( codes, K = 3 )
# back to base of "data"
setwd( '../../..' )

# now real data
setwd( 'tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01/kinship' )
message( 'Loading real data...' )
data_real <- load_pcs( codes, K = 2 )
# back to base of "data"
setwd( '../..' )


##################
### PROCESSING ###
##################

# will have to remove PC1 for popkins (and true if it were there)
# but first want to verify reason is it's highly colinear with the intercept!
n_sim <- nrow( data_sim$popkin_rom )
# construct intercept normalized as eigenvector
intercept_norm <- rep.int( 1, n_sim ) / sqrt( n_sim )
message( 'Sim: Intercept projection to PC1 of popkin ROM: ', abs( intercept_norm %*% data_sim$popkin_rom[ , 1] ) )
message( 'Sim: Intercept projection to PC1 of popkin MOR: ', abs( intercept_norm %*% data_sim$popkin_mor[ , 1] ) )

# for comparison, repeat exactly the same but on real data (where PC1 was not problematic)
n_real <- nrow( data_real$popkin_rom )
# construct intercept normalized as eigenvector
intercept_norm <- rep.int( 1, n_real ) / sqrt( n_real )
message( 'Real: Intercept projection to PC1 of popkin ROM: ', abs( intercept_norm %*% data_real$popkin_rom[ , 1] ) )
message( 'Real: Intercept projection to PC1 of popkin MOR: ', abs( intercept_norm %*% data_real$popkin_mor[ , 1] ) )


# make hacky version of data with original PC1 removed for true/popkin (shift PC2-3 down)
# this is weird for popkin MOR (PCs don't agree very well with other paired MOR bias types), but verified there's no better fit
for ( code in c('popkin_mor', 'popkin_rom') )
    data_sim[[ code ]] <- data_sim[[ code ]][ , -1 ]
# now keep only PC1-2 for all (after shifting down PCs for above methods)
data_sim <- lapply( data_sim, function( x ) x[ , 1:2 ] )

# make sure data is in the same order (as eigenvectors have arbitrary signs)
# best to align early, before calculating ranges
data_sim <- align_pcs( data_sim )

# real data only needs aligning
data_real <- align_pcs( data_real )

############
### PLOT ###
############

width <- fig_width()
fig_start(
    'pcs',
    width = width,
    height = width,
    mar_t = 1.5,
    mar_l = 4.5
)
# fill by columns!
par( mfcol = c(2,2) )

# get common ranges
xlim_sim <- range( sapply( data_sim, function( x ) range( x[ , 1 ] ) ) )
ylim_sim <- range( sapply( data_sim, function( x ) range( x[ , 2 ] ) ) )

plot_pcs( data_sim, kinship_methods$code[ kinship_methods$type == 'ROM est.' ], xlim_sim, ylim_sim, main = 'Admixed Family sim.', leg = TRUE )
mtext( 'ROM est.', side = 2, line = 3.5 )
plot_pcs( data_sim, kinship_methods$code[ kinship_methods$type == 'MOR est.' ], xlim_sim, ylim_sim )
mtext( 'MOR est.', side = 2, line = 3.5 )

# repeat for real
# get common ranges
xlim_real <- range( sapply( data_real, function( x ) range( x[ , 1 ] ) ) )
ylim_real <- range( sapply( data_real, function( x ) range( x[ , 2 ] ) ) )

plot_pcs( data_real, kinship_methods$code[ kinship_methods$type == 'ROM est.' ], xlim_real, ylim_real, main = '1000 Genomes' )
plot_pcs( data_real, kinship_methods$code[ kinship_methods$type == 'MOR est.' ], xlim_real, ylim_real )

fig_end()
