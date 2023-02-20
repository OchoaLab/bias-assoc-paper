# local links to figures
mkdir figures/
cd figures
ln -s ../../../data/sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes/kinship.pdf kinship_sim.pdf
ln -s ../../../data/sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes/auc.pdf auc_sim.pdf
ln -s ../../../data/sim-admix-n1000-m100000-k3-f0.3-s0.5-mc100-h0.8-g20-fes/pvals_eq.pdf pvals_eq_sim.pdf
ln -s ../../../data/tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01/kinship.pdf kinship_real.pdf
ln -s ../../../data/tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01/auc.pdf auc_real.pdf
ln -s ../../../data/pcs.pdf pcs.pdf
ln -s ../../../data/emin-cut.pdf emin-cut.pdf
cd ..

# separate supplement from pre-existing, preprint copy
# for this paper it's all just figures, and there's no hyperlinks worth preserving
pdftk ../bias-assoc.pdf cat 44-65 output supp.pdf
# the links remain, though they're all broken!  this removes them...
# https://askubuntu.com/questions/106154/open-source-command-line-tools-to-remove-hyperlinks-in-pdfs
sed '/Link/d' < supp.pdf > supp2.pdf
# manually compared, agreed this solves problem, so replace
mv supp2.pdf supp.pdf
# final name
mv supp.pdf File\ S1.pdf

# create zip file to upload
zip -r paper.zip paper.{pdf,tex,log} gsajnl.cls styles/ figures/ bibliography.bib genetics.bst GENETICSturquoise.pdf 
