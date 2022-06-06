### SIM ###

#m_loci=10000 # original
m_loci=100000 # bigger, matches pca-assoc
n_rep=100
gen=20

# dir output name for this run, which gets passed to other scripts
name='sim-admix-n1000-m'$m_loci'-k3-f0.3-s0.5-mc100-h0.8-g'$gen'-fes'

for rep in $(seq 1 $n_rep); do
    # simulate genotypes and phenotypes data on standard K=3 admixture
    time Rscript sim-00-sim-gen-phen.R -g $gen --fes -r $rep -m $m_loci
    # 0m10.618s small m viiiaR5
    # 0m53.628s big m

    # create all kinship estimates
    time Rscript sim-01-kinship.R --bfile $name -r $rep
    # 0m3.313s small m
    # 0m13.995s big m
    
    # run association tests
    time Rscript sim-03-assoc.R --bfile $name -p 2 -r $rep
    # 0m26.525s small m
    # 1m29.778s big m

    # calculate AUCs and SRMSDs
    time Rscript sim-04-auc-calc.R --bfile $name -r $rep
    # 0m0.615s small m
    # 0m1.487s big m
done

# kinship plots (rep-1 only)
time Rscript sim-02-kinship-plot.R --bfile $name -r 1
# 0m3.002s

# plot AUCs and SRMSDs (all reps)
time Rscript sim-05-auc-rmsd-plot.R --bfile $name --n_rep $n_rep

# statistic correlation heatmaps (all reps)
time Rscript sim-06-stats-corr.R --bfile $name --n_rep $n_rep
# 1m7.825s



### REAL ###

# 1000 Genomes version
name='tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/tgp-nygc'

cd ../data/
mkdir $name
cd $name

ln -s "$DATA_DIR/$name.bed" data.bed
ln -s "$DATA_DIR/$name.bim" data.bim
ln -s "$DATA_DIR/$name.fam" data.fam
ln -s "$DATA_DIR/pops-annot.txt" pops-annot.txt

# return to scripts dir
cd ../../scripts/

# need all kinship estimates first!
# for real data no need to do for each replicate, these are shared across replicates
# (popkin needed by simtrait)
time Rscript sim-01-kinship.R --bfile $name -r 0
# 14m6.050s/58m6.979s viiiaR5

# popkin plots! (no need to specify replicate)
time Rscript real-00-kinship-plot.R --bfile $name
# 0m9.135s

for rep in $(seq 1 $n_rep); do
    # simulate trait now!
    time Rscript real-01-simtrait.R --bfile $name --fes -r $rep
    # m_causal: 250
    # 0m2.075s

    # run association tests
    time Rscript sim-03-assoc.R --bfile $name -p 10 -r $rep
    # 108m24.091s

    # calculate AUCs and SRMSDs
    time Rscript sim-04-auc-calc.R --bfile $name -r $rep
    # 0m0.615s
done

# plot AUCs and SRMSDs (all reps)
time Rscript sim-05-auc-rmsd-plot.R --bfile $name --n_rep $n_rep
# 0m6.258s

# statistic correlation heatmaps (all reps)
time Rscript sim-06-stats-corr.R --bfile $name --n_rep $n_rep
# 8m32.654s
# NOTE: beta fig fails every time (betas are extremely correlated) 



### COMBINED ###

# compare popkin ROM/MOR (and true, in sim) in sim and TGP
time Rscript all-popkin-mor-rom-bias.R
# 0m17.224s ideapad
# remove earlier version if redoing
rm ../data/popkin-mor-rom-bias.png
# flatten huge output PDF into PNG
time pdf2png ../data/popkin-mor-rom-bias.pdf 
# 0m56.145s ideapad
# remove redundant PDF if done
rm ../data/popkin-mor-rom-bias.pdf


### THEORY/OBSOLETE ###

# empirically test the hypothesis that centering and EVD approximation commute
# this script simulates kinship as usual, we use r = K
# RESULTS:
# - commutability is not exact, especially noticeable under family structure, but still appears to be a very good approximation
# - the intercept is in the rowspace in both cases, as expected
# - rowspaces are a mess (centering does alter rowspaces enough that they really don't match; see last line of each run below, calculated via matrix ranks on concatenations of matrices)
# 
# 1) version without family structure (truly low-dimensional, K=3, Fst=0.3)
Rscript test-00-rank-pca-center.R
# RMSD center/dim-r, vs dim-r/center: 0.000501087785158289
# Corr center/dim-r, vs dim-r/center: 0.999993595141344
# Sum of kinship, center/dim-r: 1.3e-14
# Sum of kinship, dim-r/center: -3.7e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank dim-r/center + dim-r: 4
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7
#
# 2) larger K, lower FST, same result
Rscript test-00-rank-pca-center.R -k 10 -f 0.1
# RMSD center/dim-r, vs dim-r/center: 0.0005038297388324
# Corr center/dim-r, vs dim-r/center: 0.999923201852988
# Sum of kinship, center/dim-r: -1.95e-15
# Sum of kinship, dim-r/center: 1.9e-11
# Rank center/dim-r: 10
# Rank dim-r/center: 10
# Rank center/dim-r + dim-r/center: 19
# Rank center/dim-r + dim-r: 20
# Rank dim-r: 10
# Rank dim-r/center: 10
# Rank dim-r/center + dim-r: 11
# Rank intercept + dim-r: 11
# Rank intercept + dim-r/center: 11
# Rank intercept + center/dim-r: 11
# Rank intercept + dim-r/center + dim-r: 11
# Rank intercept + center/dim-r + dim-r: 20
#
# 3) back to default (K=3, Fst=0.3), but add family structure
Rscript test-00-rank-pca-center.R -g 20
# RMSD center/dim-r, vs dim-r/center: 0.00309263509709354
# Corr center/dim-r, vs dim-r/center: 0.99962525703413
# Sum of kinship, center/dim-r: -3.28e-14
# Sum of kinship, dim-r/center: 4.45e-11
# Rank center/dim-r: 3
# Rank dim-r/center: 3
# Rank center/dim-r + dim-r/center: 6
# Rank center/dim-r + dim-r: 6
# Rank dim-r: 3
# Rank dim-r/center: 3
# Rank intercept + dim-r: 4
# Rank intercept + dim-r/center: 4
# Rank intercept + center/dim-r: 4
# Rank intercept + dim-r/center + dim-r: 4
# Rank intercept + center/dim-r + dim-r: 7

