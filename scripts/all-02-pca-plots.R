
library(genio)       # to write BED files for external software



## real dataset (1000 Genomes)
setwd( 'D:/3.Duke/research/alex_ochoa/1215_2021_PCA/kinship-real' )


## read data
kinship_gcta <- read_grm( 'gcta' )$kinship

kinship_popkin <- read_grm( 'popkin' )$kinship

kinship_std_mor <- read_grm( 'std_mor' )$kinship

kinship_std_rom <- read_grm( 'std_rom' )$kinship

kinship_wg <- read_grm( 'wg' )$kinship



#Calculate PCs using the eigendecomposition of the kinship matrix:

pcs_gcta <- eigen( kinship_gcta )$vectors

pcs_popkin <- eigen( kinship_popkin )$vectors

pcs_std_mor <- eigen( kinship_std_mor )$vectors

pcs_std_rom <- eigen( kinship_std_rom )$vectors

pcs_wg <- eigen( kinship_wg )$vectors



# PCs are the columns, so to plot the first two PCs:

plot( pcs_gcta[,1], pcs_gcta[,2],xlab='PC1',ylab='PC2',main='gcta-real' , xlim=c(-0.03,0.04), ylim=c(-0.03,0.04))

plot( pcs_popkin[,1], pcs_popkin[,2],xlab='PC1',ylab='PC2',main='popkin-real', xlim=c(-0.03,0.04), ylim=c(-0.03,0.04) )

plot( pcs_std_mor[,1], pcs_std_mor[,2],xlab='PC1',ylab='PC2',main='std_mor-real', xlim=c(-0.03,0.04), ylim=c(-0.03,0.04) )

plot( pcs_std_rom[,1], pcs_std_rom[,2],xlab='PC1',ylab='PC2',main='std_rom-real', xlim=c(-0.03,0.04), ylim=c(-0.03,0.04) )

plot( pcs_wg[,1], pcs_wg[,2],xlab='PC1',ylab='PC2',main='wg-real', xlim=c(-0.03,0.04), ylim=c(-0.03,0.04) )

## all on one

plot( pcs_gcta[,1], pcs_gcta[,2],xlab='PC1',ylab='PC2' , xlim=c(-0.03,0.04), ylim=c(-0.03,0.04),pch = 1, col = 1)

points( pcs_popkin[,1], pcs_popkin[,2],pch = 2, col = 2) 

points( pcs_std_mor[,1], pcs_std_mor[,2],pch = 3, col = 3)

points( pcs_std_rom[,1], pcs_std_rom[,2],pch = 4, col = 4)

points( pcs_wg[,1], pcs_wg[,2],pch = 5, col = 5)

legend(0.02,0.045,c("gcta","popkin","std_mor","std_rom","WG"),pch=1:5,col=1:5,bty = "n")

## admixed family simulation

setwd( 'D:/3.Duke/research/alex_ochoa/1215_2021_PCA/kinship-sim' )


## read data
kinship_gcta <- read_grm( 'gcta' )$kinship

kinship_gcta_lim <- read_grm( 'gcta_lim' )$kinship

kinship_popkin <- read_grm( 'popkin' )$kinship

kinship_true <- read_grm( 'true' )$kinship

kinship_std_mor <- read_grm( 'std_mor' )$kinship

kinship_std_rom <- read_grm( 'std_rom' )$kinship

kinship_std_rom_lim <- read_grm( 'std_rom_lim' )$kinship

kinship_wg <- read_grm( 'wg' )$kinship

kinship_wg_lim <- read_grm( 'wg_lim' )$kinship


#Calculate PCs using the eigendecomposition of the kinship matrix:

pcs_gcta <- eigen( kinship_gcta )$vectors
pcs_gcta_lim <- eigen( kinship_gcta_lim )$vectors

pcs_popkin <- eigen( kinship_popkin )$vectors
pcs_true <- eigen( kinship_true )$vectors

pcs_std_mor <- eigen( kinship_std_mor )$vectors
pcs_std_rom <- eigen( kinship_std_rom )$vectors
pcs_std_rom_lim <- eigen( kinship_std_rom_lim )$vectors

pcs_wg <- eigen( kinship_wg )$vectors
pcs_wg_lim <- eigen( kinship_wg_lim )$vectors



# PCs are the columns, so to plot the first two PCs:

plot( pcs_gcta[,1], pcs_gcta[,2],xlab='PC1',ylab='PC2',main='gcta-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))
plot( pcs_gcta_lim[,1], pcs_gcta_lim[,2],xlab='PC1',ylab='PC2',main='gcta-lim-sim',xlim=c(-0.05,0.065), ylim=c(-0.06,0.065) )


plot( pcs_popkin[,1], pcs_popkin[,2],xlab='PC1',ylab='PC2',main='popkin-sim',xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))
plot( pcs_true[,1], pcs_true[,2],xlab='PC1',ylab='PC2',main='true-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065) )

plot( pcs_std_mor[,1], pcs_std_mor[,2],xlab='PC1',ylab='PC2',main='std_mor-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))
plot( pcs_std_rom[,1], pcs_std_rom[,2],xlab='PC1',ylab='PC2',main='std_rom-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))
plot( pcs_std_rom_lim[,1], pcs_std_rom_lim[,2],xlab='PC1',ylab='PC2',main='std_rom_lim-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))

plot( pcs_wg[,1], pcs_wg[,2],xlab='PC1',ylab='PC2',main='wg-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))
plot( pcs_wg_lim[,1], pcs_wg_lim[,2],xlab='PC1',ylab='PC2',main='wg_lim-sim' ,xlim=c(-0.05,0.065), ylim=c(-0.06,0.065))

## all on one

plot( pcs_true[,1], pcs_true[,2],xlab='PC1',ylab='PC2' , xlim=c(-0.05,0.065), ylim=c(-0.06,0.07),pch = 1, col = 1)

points( pcs_gcta[,1], pcs_gcta[,2],pch = 2, col = 2) 

points( pcs_gcta_lim[,1], pcs_gcta_lim[,2],pch = 3, col = 3)

points( pcs_popkin[,1], pcs_popkin[,2],pch = 4, col = 4)

points( pcs_std_mor[,1], pcs_std_mor[,2],pch = 5, col = 5)

points( pcs_std_rom[,1], pcs_std_rom[,2],pch = 6, col = 6)

points( pcs_std_rom_lim[,1], pcs_std_rom_lim[,2],pch = 7, col = 7)

points( pcs_wg[,1], pcs_wg[,2],pch = 8, col = 8)

points( pcs_wg_lim[,1], pcs_wg_lim[,2],pch = 9, col =11)


legend(0.045,0.08,c("true","gcta","gcta_lim","popkin","std_mor","std_rom","std_rom_lim","WG","wg_lim"),
       pch=1:9,col=c(1:8,11),bty = "n")

