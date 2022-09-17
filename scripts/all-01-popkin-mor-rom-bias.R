library(genio) # read_grm
library(popkin) # inbr
library(ochoalabtools) # fig_start

# go where the data is
setwd( '../data/' )

# start figure in this location
scale <- fig_scale( 2 ) # w/h ratio
fig_start(
    'popkin-mor-rom-bias',
    width = scale[1],
    height = scale[2],
    mar_t = 1.5
)
# two panels
par( mfrow = c(1, 2) )

# simulation first
setwd( 'sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes' )

# use first replicate only
setwd( 'rep-1/kinship' )

# load the three unbiased-type kinship matrices
# transform diagonals to inbreeding immediately
kinship_true <- inbr_diag( read_grm( 'true' )$kinship / 2 )
kinship_popkin_rom <- inbr_diag( read_grm( 'popkin_rom' )$kinship / 2 )
kinship_popkin_mor <- inbr_diag( read_grm( 'popkin_mor' )$kinship / 2 )
# set range of plot, symmetric so its more pleasing
kin_max <- max( kinship_true, kinship_popkin_rom, kinship_popkin_mor )

# initialize plot area
plot(
    NA,
    xlim = c(0, kin_max),
    ylim = c(0, kin_max),
    xlab = 'True kinship',
    ylab = 'Popkin estimates',
    main = 'Admixed Family sim.'
)
panel_letter( 'A', adj = -0.2 )
# plot y=x diagonal in background
abline( 0, 1, lty = 2, col = 'grey' )

# plot off-diagonal without duplicating dots
indexes <- lower.tri( kinship_true )
points( kinship_true[ indexes ], kinship_popkin_rom[ indexes ], pch = '.' )
# diagonal use diff color
points( diag( kinship_true ), diag( kinship_popkin_rom ), pch = '.', col = 'blue' )

# repeat MOR (will be on top of ROM)
points( kinship_true[ indexes ], kinship_popkin_mor[ indexes ], pch = '.', col = 'red' )
points( diag( kinship_true ), diag( kinship_popkin_mor ), pch = '.', col = 'orange' )

legend(
    'topleft',
    c( 'ROM kinship', 'ROM inbreeding', 'MOR kinship', 'MOR inbreeding' ),
    col = c( 'black', 'blue', 'red', 'orange' ),
    pch = 19,
    bty = 'n',
    cex = 0.7
)


# now real data
setwd( '../../../tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01/kinship' )

# load the two unbiased-type kinship matrix estimates (no "true" for real data)
kinship_popkin_rom <- inbr_diag( read_grm( 'popkin_rom' )$kinship / 2 )
kinship_popkin_mor <- inbr_diag( read_grm( 'popkin_mor' )$kinship / 2 )

# set range of plot, symmetric so its more pleasing
kin_max <- max( kinship_popkin_rom, kinship_popkin_mor )

# initialize plot area
plot(
    NA,
    xlim = c(0, kin_max),
    ylim = c(0, kin_max),
    xlab = 'Popkin ROM est.',
    ylab = 'Popkin MOR est.',
    main = '1000 Genomes'
)
panel_letter( 'B', adj = -0.2 )
# plot y=x diagonal in background
abline( 0, 1, lty = 2, col = 'grey' )

# plot off-diagonal without duplicating dots
indexes <- lower.tri( kinship_popkin_rom )
points( kinship_popkin_rom[ indexes ], kinship_popkin_mor[ indexes ], pch = '.', col = 'red' )
# diagonal use diff color?
points( diag( kinship_popkin_rom ), diag( kinship_popkin_mor ), pch = '.', col = 'orange' )

legend(
    'topleft',
    c( 'Kinship', 'Inbreeding' ),
    col = c( 'red', 'orange' ),
    pch = 19,
    bty = 'n',
    cex = 0.7
)


fig_end()

