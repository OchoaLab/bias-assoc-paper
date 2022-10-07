# run LMM assoc in R, with detailed reporting to understand what is happening in reality for non-posdef kinship
# will focus on true vs lims for simplicity

library(optparse)
library(readr)
library(ochoalabtools)

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

# go where we want data outputs to be
setwd( '../data' )
setwd( dir_out )
dir_rep <- paste0( 'rep-', rep )
setwd( dir_rep )

# load precalculated data inside rep
# the data we calculated in R using GLS with V's from REML
data <- read_tsv( 'lmm-intercept-test.txt.gz', col_types = 'icdd' )
# the data directly and entirely calculated by GCTA
betas <- read_tsv( 'betas.txt.gz', show_col_types = FALSE )

# start by validating our calculations against GCTA's betas (they don't provide intercept)
codes <- unique( data$kinship )
for ( code in codes ) {
    # get GCTA's betas for this matrix
    code_betas <- paste0( 'lmm_', code )
    betas_gcta <- betas[[ code_betas ]]
    # and our betas
    betas_gls <- data$beta[ data$kinship == code ]
    # visually they're a great match!
    # make sure mean abs error is small enough
    stopifnot( mean( abs( betas_gcta - betas_gls ) ) < 1e-4 )
}
# cool, we're all set in terms of agreeing with GCTA!  No need for that data anymore

# now compare our betas internally, which visually also agree to a great extent
# because we didn't truncate digits as much as GCTA, here agreements are even greater!
for ( i in 2L : length( codes ) ) {
    code_i <- codes[ i ]
    for ( j in 1L : ( i - 1L ) ) {
        code_j <- codes[ j ]
        diff_ij <- mean( abs( data$beta[ data$kinship == code_i ] - data$beta[ data$kinship == code_j ] ) )
        #message( code_i, ', ', code_j, ': ', diff_ij )
        stopifnot( diff_ij < 1e-7 )
    }
}

# excellent, so previous theory and results (concerning betas) is largely validated here!
# now the intercepts which are the new question we haven't been able to probe before empirically
alpha_tr <- data$intercept[ data$kinship == 'true' ]
alpha_wg <- data$intercept[ data$kinship == 'wg_rom_lim' ]
alpha_st <- data$intercept[ data$kinship == 'std_rom_lim' ]

# for intercepts there are two cases:
# true and WG appear to agree perfectly!  (supported by theory)
# confirmed visually and here we'll just test as a precision test
stopifnot( mean( abs( alpha_tr - alpha_wg ) ) < 1e-7 )

# save fig in lower level
setwd( '..' )

# STD actually differs from true's intercept, as expected from original theory
# the plot is very strange, with a sort of funnel in the middle where agreement is greater
width <- fig_width() / 2
fig_start(
    'lmm-intercept-test',
    width = width,
    height = width
)
plot(
    alpha_tr,
    alpha_st - alpha_tr,
    xlab = expression( bold( alpha ) ),
    ylab = expression( bold( eta == alpha*minute - alpha ) ),
    pch = '.'
)
fig_end()

## plot(
##     alpha_tr,
##     alpha_wg - alpha_tr
## )
## abline( h = 0, lty = 2, col = 'gray' )
