library(optparse) # for terminal options
library(readr)    # to read tables
library(popkin)   # to plot
library(ochoalabtools) # for nice PDF

# constants
tolp <- 1e-2
tolb <- 1e-3

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Directory to process (under ../data/, containing input plink files data.BED/BIM/FAM/PHEN)", metavar = "character"),
    make_option("--n_rep", type = "integer", default = NA, 
                help = "Total number of replicates", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
dir_out <- opt$bfile
n_rep <- opt$n_rep
if ( is.na( n_rep ) )
    stop( 'Option `--n_rep` is required!' )

# before switching away from "scripts", load a table located there
kinship_methods <- read_tsv( 'kinship_methods.txt', col_types = 'cccc' )

# load pre-existing data
setwd( '../data/' )
setwd( dir_out )

# if this directory exists here, this is real data
is_sim <- !file.exists( 'kinship' )

# in real data the oracle methods (truth and biased limits) are missing, exclude from `kinship_methods`
if ( !is_sim )
    kinship_methods <- kinship_methods[ kinship_methods$type != 'ROM lim.', ]

# let's plot data in the order of the `kinship_methods` table
# also, plot PCA first, LMM second
method_codes <- c(
    paste0( 'pca_', kinship_methods$code ),
    paste0( 'lmm_', kinship_methods$code )
)

# compute sum of equal elements for every pair of columns (like a correlation), but up to a given tolerance
pair_identical <- function( X, tol ) {
    n <- ncol( X )
    Y <- matrix( 0, n, n )
    # can't think of a better way than actually navigating all pairs
    for ( i in 2 : n ) {
        xi <- X[ , i ]
        # diagonal that makes sense is number of self non-misisng elements
        Y[ i, i ] <- sum( !is.na( xi ) )
        for ( j in 1 : (i-1) ) {
            xj <- X[ , j ]
            # this is the calculation we were interested in
            Y[ i, j ] <- sum( abs( xi - xj ) < tol, na.rm = TRUE )
            # copy both ways
            Y[ j, i ] <- Y[ i, j ]
        }
    }
    # above loop skips i=1, do now
    Y[ 1, 1 ] <- sum( !is.na( X[,1] ) )
    return( Y )
}


process_tables <- function( name, tol ) {
    # read files
    file <- paste0( name, '.txt.gz' )
    data <- read_tsv( file, col_types = cols( ) )
    # reorder columns of data to be a manually-selected order
    data <- data[ method_codes ]
    # rest makes most sense as matrix
    data <- as.matrix( data )
    # compute correlation matrix
    C <- cor( data, use = 'pairwise.complete.obs' )

    # replicate calculation but with the equality comparisons
    # denominator first, sum of non-missing cases
    M <- crossprod( !is.na( data ) )
    # numerator is more complicated...
    E <- pair_identical( data, tol = tol )

    # return all three parts
    return( list(
        C = C,
        M = M,
        E = E
    ) )
}

# get data from each replicate
# add to various running sums
for ( rep in 1 : n_rep ) {
    setwd( paste0( 'rep-', rep ) )

    message( 'rep-', rep, ': pvals' )
    data <- process_tables( 'pvals', tol = tolp )
    if ( rep == 1 ) {
        Cp <- data$C
        Mp <- data$M
        Ep <- data$E
    } else {
        Cp <- Cp + data$C
        Mp <- Mp + data$M
        Ep <- Ep + data$E
    }

    # repeat for betas
    message( 'rep-', rep, ': betas' )
    data <- process_tables( 'betas', tol = tolb )
    if ( rep == 1 ) {
        Cb <- data$C
        Mb <- data$M
        Eb <- data$E
    } else {
        Cb <- Cb + data$C
        Mb <- Mb + data$M
        Eb <- Eb + data$E
    }

    setwd( '..' )
}

# complete averages
# for correlations it's a naive average (wrong considering missingess, but meh)
Cp <- Cp / n_rep
Cb <- Cb / n_rep
# for equality, properly consider missingness
Ep <- Ep / Mp
Eb <- Eb / Mb

# save all of these tables!
write_tsv( as.data.frame( Cp ), 'pvals_cor.txt' )
write_tsv( as.data.frame( Cb ), 'betas_cor.txt' )
write_tsv( as.data.frame( Ep ), 'pvals_eq.txt' )
write_tsv( as.data.frame( Eb ), 'betas_eq.txt' )
