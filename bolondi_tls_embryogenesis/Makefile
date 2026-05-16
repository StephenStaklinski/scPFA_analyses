export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
PYPLOTTING_SIF := ${CONTAINERS}/pyplotting/pyplotting.sif
SCANPY_SIF := ${CONTAINERS}/scanpy/scanpy.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
MONOCLE_SIF := ${CONTAINERS}/monocle/monocle.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
GEXLINEAGE := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexLineage

CURR_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis

DATA_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis/data/GSE220949
TLS1_PATH := $(DATA_DIR)/GSM6831757_TracerTLS1
TLS2_PATH := $(DATA_DIR)/GSM6831759_TracerTLS2
TLSMULTI_PATH := $(DATA_DIR)/GSM6831760_TracerTLSMulti

CHARMATRICES := tls1_character_matrix.tsv tls2_character_matrix.tsv
GEXMATS := tls1_gex.tsv tls2_gex.tsv
VINETREES := tls1.vine.nex tls2.vine.nex
GEXLINEAGERESULTS := tls1.fit.summary.tsv tls2.fit.summary.tsv
ORIGEMBEDDINGS := tls1_gex.umap.pdf tls2_gex.umap.pdf
PALANTIRTRAJECTORIES := tls1_gex.palantir.pseudotime.tsv tls2_gex.palantir.pseudotime.tsv
MONOCLETRAJECTORIES := tls1_gex.monocle.pseudotime.pdf tls2_gex.monocle.pseudotime.pdf
LOGPLOTS := tls1.fit.log.pdf tls2.fit.log.pdf
PANELPLOTS := tls1.fit.panel.png tls2.fit.panel.png
GEXLINEAGEGSEA := tls1.fit.L.gsea.tsv tls2.fit.L.gsea.tsv
FACTORHEATMAPS := tls1.fit.L.factor_heatmap.pdf tls2.fit.L.factor_heatmap.pdf
FLOWTREES := tls1.fit.latent_flow.factor_trees.circular.pdf tls2.fit.latent_flow.factor_trees.circular.pdf
EUCLIDEANPSEUDOTIME := tls1.fit.euclidean_factor_pseudotime.pdf tls2.fit.euclidean_factor_pseudotime.pdf

all: prep_data run_vine

prep_data: $(CHARMATRICES) $(GEXMATS) $(ORIGEMBEDDINGS)
palantir: $(PALANTIRTRAJECTORIES)
monocle: $(MONOCLETRAJECTORIES)
run_vine: $(VINETREES)
gexlineage: $(GEXLINEAGERESULTS)
plotting: $(LOGPLOTS) $(PANELPLOTS) $(GEXLINEAGEGSEA) $(FACTORHEATMAPS) $(FLOWTREES) $(EUCLIDEANPSEUDOTIME)

# Obtain the TLS1 character matrix
tls1_character_matrix.tsv tls1_mutation_dict.tsv: $(TLS1_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls1_character_matrix.tsv /mnt/cwd/tls1_mutation_dict.tsv > tls1_character_matrix.term

# Obtain the TLS2 character matrix
tls2_character_matrix.tsv tls2_mutation_dict.tsv: $(TLS2_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls2_character_matrix.tsv /mnt/cwd/tls2_mutation_dict.tsv > tls2_character_matrix.term

# TODO: Split the TLSMulti allele table into separate tables for each of the additional 21 samples, obtain the character matrices, and fill in the other steps below

# Subset the TLS1 gene expression data to the cells present in the character matrix
tls1_gex.tsv: $(TLS1_PATH)/matrix.mtx $(TLS1_PATH)/barcodes.tsv $(TLS1_PATH)/features.tsv tls1_character_matrix.tsv
		singularity exec --bind $(TLS1_PATH):/mnt/tls1 --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/tls1_character_matrix.tsv /mnt/tls1/matrix.mtx \
			/mnt/tls1/barcodes.tsv /mnt/tls1/features.tsv /mnt/cwd/tls1_gex.tsv > tls1_gex.term

# Subset the TLS2 gene expression data to the cells present in the character matrix
tls2_gex.tsv: $(TLS2_PATH)/matrix.mtx $(TLS2_PATH)/barcodes.tsv $(TLS2_PATH)/features.tsv tls2_character_matrix.tsv
		singularity exec --bind $(TLS2_PATH):/mnt/tls2 --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/tls2_character_matrix.tsv /mnt/tls2/matrix.mtx \
			/mnt/tls2/barcodes.tsv /mnt/tls2/features.tsv /mnt/cwd/tls2_gex.tsv > tls2_gex.term

# Plot the original gex embeddings
%_gex.pca.pdf %_gex.umap.pdf %_gex.pca.tsv %_gex.umap.tsv: %_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(SCANPY_SIF) \
		python /mnt/src/scanpy_dim_reduction.py /mnt/cwd/$*_gex.tsv /mnt/cwd/$*_gex

# Get gex trajectories with palantir
%_gex.palantir.pseudotime.tsv: %_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(SCANPY_SIF) \
		python /mnt/src/palantir_trajectory.py /mnt/cwd/$*_gex.tsv /mnt/cwd/$*_gex

# Get gex trajectories with monocle
%_gex.monocle.pseudotime.pdf: %_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(MONOCLE_SIF) \
		Rscript /mnt/src/monocle_trajectory.R /mnt/cwd/$*_gex.tsv /mnt/cwd/$*_gex

# Run VINE on the character matrices to infer lineage trees
%.vine.nwk: %_character_matrix.tsv
	/usr/bin/time -o $*.vine.time $(VINE) -M 400 -c 100 -v 0 \
	-i CRISPR $*_character_matrix.tsv --logf $*.vine.log --mean $*.vine.mean.nwk > $*.vine.nwk

# Convert nwk to nexus
%.vine.nex: %.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd/ --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py /mnt/cwd/$*.vine.nwk /mnt/cwd/$*.vine.nex

# Fit the gene expression model
%.fit.summary.tsv %.fit.F.tsv %.fit.L.tsv %.fit.X.tsv %.fit.log %.fit.latent_flow.tsv: %.vine.nex %_gex.tsv
	/usr/bin/time -o $*.fit.time ${GEXLINEAGE} \
	--trees $*.vine.nex \
	--expr $*_gex.tsv \
	--tree-total-time 1.0 \
	--outprefix $*.fit \
	--filter-test moran \
	--dim 10 \
	--varimax \
	--L-l1-strength 1.0 \
	--max-iter 10000 \
	--nthreads 1 \
	--verbose-log \
	--seed 1 > $*.fit.term

# Plot the fit optimization log results
%.fit.log.pdf: %.fit.log
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_optimization_log.py /mnt/cwd/$*.fit.log -1 0 /mnt/cwd/$*.fit.log

# Plot the fit results
%.fit.panel.png: %.fit.F.tsv %.fit.L.tsv %.fit.X.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_fit_panel.py \
			/mnt/cwd/$*.fit.F.tsv \
			/mnt/cwd/$*.fit.L.tsv \
			/mnt/cwd/$*.fit.X.tsv \
			/mnt/cwd/$*.fit.panel.png

# Run GO BP GSEA on fitted L matrix
%.fit.L.gsea.tsv: %.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) \
		/mnt/src/run_factor_gsea.R /mnt/cwd/$*.fit.L.tsv mouse /mnt/cwd/$*.fit.L.gsea > $*.fit.gsea.term

# Plot the factor heatmap
%.fit.L.factor_heatmap.pdf: %.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_heatmap.py \
			/mnt/cwd/$*.fit.L.tsv \
			/mnt/cwd/$*.fit.L.factor_heatmap.pdf

# Plot the latent flow results
%.fit.latent_flow.factor_trees.pdf %.fit.latent_flow.total_euclidean_distance.tsv: %.fit.latent_flow.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow.py \
			/mnt/cwd/$*.fit.latent_flow.tsv \
			/mnt/cwd/$*.fit.latent_flow

# Plot euclidean latent-factor pseudotime on gex-only PCA/UMAP embeddings
%.fit.euclidean_factor_pseudotime.pdf: %.fit.latent_flow.total_euclidean_distance.tsv %_gex.pca.tsv %_gex.umap.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_euclidean_distance_pseudotime.py \
			/mnt/cwd/$*_gex.pca.tsv \
			/mnt/cwd/$*_gex.umap.tsv \
			/mnt/cwd/$*.fit.latent_flow.total_euclidean_distance.tsv \
			/mnt/cwd/$*.fit.latent_flow

