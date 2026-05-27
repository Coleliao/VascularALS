# recording command line running pyscenic
# Kuo 2026/03/11

conda activate pyscenic
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH



tf='~/02.ref_genome/cisTarget/human/allTFs_hg38.txt'
feather='~/02.ref_genome/cisTarget/human/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather' # here use a range of +-10k bp 
tbl='~/02.ref_genome/cisTarget/human/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'

loom='./capillary_endo_counts.loom'

## step 1 GRN construction 

pyscenic grn \
	--num_workers 10 \
	--output ./results/grn_output.tsv \
	--method grnboost2 \
	$loom $tf

## step 2 TF-cis-targets to get Regulon (a TF and its all target genes)
pyscenic ctx ./results/grn_output.tsv \
	$feather \
	--annotations_fname $tbl \
	--expression_mtx_fname $loom \
	--mode 'dask_multiprocessing' \
	--output ./results/ctx_output.csv \
	--num_workers 10 \
	--mask_dropouts

# step 3 AUCell
pyscenic aucell \
	$loom \
	./results/ctx_output.csv \
	--output ./results/aucell_output.csv ./results/aucell_output.loom \
	--num_workers 10


