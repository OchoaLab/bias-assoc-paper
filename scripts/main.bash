### SIM ###

m_loci=100000
n_rep=100
gen=20

# dir output name for this run, which gets passed to other scripts
name='sim-admix-n1000-m'$m_loci'-k3-f0.3-s0.5-mc100-h0.8-g'$gen'-fes'

for rep in $(seq 1 $n_rep); do
    # simulate genotypes and phenotypes data on standard K=3 admixture
    time Rscript sim-00-sim-gen-phen.R -g $gen --fes -r $rep -m $m_loci
    # 0m53.628s viiiaR5

    # create all kinship estimates
    time Rscript sim-01-kinship.R --bfile $name -r $rep
    # 0m13.995s
    
    # run association tests
    time Rscript sim-03-assoc.R --bfile $name -p 2 -r $rep
    # 1m29.778s

    # calculate AUCs and SRMSDs
    time Rscript sim-04-auc-calc.R --bfile $name -r $rep
    # 0m1.487s
done

# kinship plots (rep-1 only)
time Rscript sim-02-kinship-plot.R --bfile $name -r 1
# 0m3.002s

# gather AUC/RMSD data from all reps into single table
time Rscript sim-05-auc-rmsd-table.R --bfile $name --n_rep $n_rep
# 0m1.922s viiiaR5

# plot AUCs and SRMSDs (all reps)
time Rscript sim-05-auc-rmsd-plot.R --bfile $name --n_rep $n_rep

# statistic correlation heatmaps (all reps)
time Rscript sim-06-stats-corr.R --bfile $name --n_rep $n_rep
# 1m7.825s

# obtain variance component estimates for many follow-up analyses
time Rscript sim-07-reml.R --bfile $name --n_rep $n_rep
# 12m23.149s/44m13.530s

# calculate V matrices for all cases, save them
time Rscript sim-08-calc-Vs.R --bfile $name --n_rep $n_rep
# 3m37.934s labbyDuke

# calculate minimum eigenvalues and condition numbers for all kinship and V matrices
time Rscript sim-09-eigen.R --bfile $name --n_rep $n_rep
# 692m24.452s/7450m57.683s viiiaR5

# generate refitted intercepts to test what happens to them
time Rscript sim-10-lmm-intercept-test.R --bfile $name -r 1
# 54m34.021s/136m14.469s ideapad
# plot non-trivial results!
time Rscript sim-11-lmm-intercept-test-plot.R --bfile $name -r 1
# creates: lmm-intercept-test.pdf
# 0m1.688s viiiaR5


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
# a version without WG better for presentations
time Rscript real-00-kinship-plot.R --bfile $name --noWG

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

# # DCC version
# # after making traits locally (fast), run association tests on cluster
# # first transfer data:
# cd ../data/$name/
# tar -chzf data.tgz data.{bed,bim,fam} pops-annot.txt kinship/ rep-*/{data.phen,simtrait.RData}
# # transfer to DCC after creating output directory structure on DCC
# scp data.tgz $dcc:/work/ao128/bias-assoc-paper/data/$name/
# rm data.tgz # cleanup

# # on DCC
# cd /work/ao128/bias-assoc-paper/data/$name/
# tar -xzf data.tgz
# rm data.tgz
# # when jobs are done, from this same location:
# tar -czf pvals.tgz rep-*/{pvals,betas}.txt.gz
# rm pvals.tgz # after transfer

# # back on computer
# cd ~/docs/ochoalab/bias-assoc-paper/data/$name/
# scp $dcc:/work/ao128/bias-assoc-paper/data/$name/pvals.tgz
# tar -xzf pvals.tgz
# rm pvals.tgz

# gather AUC/RMSD data from all reps into single table
time Rscript sim-05-auc-rmsd-table.R --bfile $name --n_rep $n_rep
# 0m2.410s viiiaR5

# plot AUCs and SRMSDs (all reps)
time Rscript sim-05-auc-rmsd-plot.R --bfile $name --n_rep $n_rep
# 0m6.258s
# a version without WG better for presentations
time Rscript sim-05-auc-rmsd-plot.R --bfile $name --n_rep $n_rep --noWG

# a hacky smaller version for grants only
time Rscript sim-05-auc-plot_lmm-small.R --bfile $name --n_rep $n_rep

# statistic correlation heatmaps (all reps)
# first calculate, save tables
time Rscript sim-06-stats-corr-calc.R --bfile $name --n_rep $n_rep
# 6m12.550s
# now plot!
time Rscript sim-06-stats-corr-plot.R --bfile $name
# 0m0.896s
# NOTE: beta fig fails every time (betas are extremely correlated) 
# a version without WG better for presentations
time Rscript sim-06-stats-corr-plot.R --bfile $name --noWG
# 0m0.863s

# obtain variance component estimates for many follow-up analyses
time Rscript sim-07-reml.R --bfile $name --n_rep $n_rep
# 308m42.172s/474m56.338s

# calculate V matrices for all cases, save them
time Rscript sim-08-calc-Vs.R --bfile $name --n_rep $n_rep
# 8m40.166s labbyDuke

# calculate minimum eigenvalues and condition numbers for all kinship and V matrices
time Rscript sim-09-eigen.R --bfile $name --n_rep $n_rep
# 13m27.778s/77m36.703s viiiaR5



### COMBINED ###

# compare popkin ROM/MOR (and true, in sim) in sim and TGP
time Rscript all-01-popkin-mor-rom-bias.R
# 0m17.224s ideapad
# remove earlier version if redoing
rm ../data/popkin-mor-rom-bias.png
# flatten huge output PDF into PNG
time pdf2png ../data/popkin-mor-rom-bias.pdf 
# 0m56.145s ideapad
# remove redundant PDF if done
rm ../data/popkin-mor-rom-bias.pdf

# create plot that compares PCs across bias types
time Rscript all-02-pca-plots.R 
# 0m21.909s viiiaR5
# Sim: Intercept projection to PC1 of popkin ROM: 0.939904048930533
# Sim: Intercept projection to PC1 of popkin MOR: 0.896792315879261
# Real: Intercept projection to PC1 of popkin ROM: 0.89979532627367
# Real: Intercept projection to PC1 of popkin MOR: 0.955378629802354

# plots variance component values and prediction errors
time Rscript all-03-reml-pred-sigmas.R
# 0m0.886s labbyDuke

# plot minimum eigenvalues, non-posdef proportions, and condition number distributions
time Rscript all-04-emin-kappa.R
# creates: emin.pdf, emin-cut.pdf, kappa.pdf
# 0m0.670s viiiaR5

# makes figures that proves several things:
# - for admix sim
#   - RMSD and AUC deviations form expectation (WG) are perfectly predicted by kappa
#   - sigma errors are uncorrelated to kappa, so errors are restricted to assoc step
# - for TGP
#   - RMSD and AUC errors are generally smaller but still extrema are not predicted at all by kappa
#   - sigma errors are way biggger and also not predicted by kappa
#   - sigma errors drive outliers here and not assoc step
time Rscript all-05-kappa-vs-pred-err.R
# creates: kappa-vs-pred-err.pdf, reml-err-vs-pred-err.pdf
# 0m0.907s viiiaR5



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

