# usage: export OG='HER7'; ../scripts/visualize_trees.sh
#OG='HER7'
SCRIPTS='../scripts'
SPECIESTREE="$SCRIPTS/SpeciesTree_model.tree"

module load conda-env/env.genomics
module load notung

echo python scripts/mpr_tree.py $OG.fast.tree 
python $SCRIPTS/mpr_tree.py $OG.fast.tree
python $SCRIPTS/visualize_tree.py $OG.fast.mpr.tree $OG.aln $SCRIPTS/species_ids.txt $SCRIPTS/lineage_colors.dmp 0.95 $OG.fast.mpr.pdf

PRENUM=`grep ">" $OG.aln | wc -l`
NEWNUM=`grep ">" $OG.trim.aln | wc -l`
PRELEN=`python $SCRIPTS/aln_len.py $OG.aln`
NEWLEN=`python $SCRIPTS/aln_len.py $OG.trim.aln`

echo -e "\nFLAG: $GENE alignment length pre trimming = $PRELEN" 
echo -e "FLAG: $GENE alignment length post trimming = $NEWLEN" 
echo -e "FLAG: $GENE seq count pre trimming = $PRENUM" 
echo -e "FLAG: $GENE seq count post trimming = $NEWNUM\n" 

echo python scripts/mpr_tree.py $OG.contree 
python $SCRIPTS/mpr_tree.py $OG.contree

echo java -jar /depot/jwisecav/apps/bell/notung/Notung-2.9.1.5.jar --root --treeoutput newick --nolosses -g $OG.mpr.tree -s $SPECIESTREE --infertransfers false --speciestag postfix --usegenedir 
java -jar /depot/jwisecav/apps/bell/notung/Notung-2.9.1.5.jar --root --treeoutput newick --nolosses -g $OG.mpr.tree -s $SPECIESTREE --infertransfers false --speciestag postfix --usegenedir

mv $OG.mpr.tree.rooting.0 $OG.ntr.tree

python $SCRIPTS/visualize_tree.py $OG.mpr.tree $OG.aln $SCRIPTS/species_ids.txt $SCRIPTS/lineage_colors.dmp 95 $OG.mpr.pdf
python $SCRIPTS/visualize_tree.py $OG.ntr.tree $OG.aln $SCRIPTS/species_ids.txt $SCRIPTS/lineage_colors.dmp 95 $OG.ntr.pdf

