export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXPCA := $(GEX_LINEAGE_DIR)/bin/gexPca
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
HOTSPOT_SIF := $(CONTAINERS)/hotspot/hotspot.sif
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_pca
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_filter
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data

CLONES := 1 2 3 4
PCA_DIM := 50
PCA_FLOW_DIMS ?= 10
PCA_EIGENVECTORS := $(CLONES:%=quinn_%.pca.eigenvectors.tsv)
PCA_TIMES := $(CLONES:%=quinn_%.pca.time)
PCA_GSEA := $(CLONES:%=quinn_%.pca.gsea.tsv)
NO_RIBO_MITO_PCA_EIGENVECTORS := $(CLONES:%=quinn_%_no_ribo_mito.pca.eigenvectors.tsv)
NO_RIBO_MITO_PCA_TIMES := $(CLONES:%=quinn_%_no_ribo_mito.pca.time)
NO_RIBO_MITO_PCA_GSEA := $(CLONES:%=quinn_%_no_ribo_mito.pca.gsea.tsv)
PCA_FLOW_STAMPS := $(CLONES:%=quinn_%.pca.latent_flow.plots.done)
NO_RIBO_MITO_PCA_FLOW_STAMPS := $(CLONES:%=quinn_%_no_ribo_mito.pca.latent_flow.plots.done)

.PHONY: all clean
all: $(PCA_EIGENVECTORS) $(PCA_TIMES) $(PCA_GSEA) \
	$(NO_RIBO_MITO_PCA_EIGENVECTORS) $(NO_RIBO_MITO_PCA_TIMES) \
	$(NO_RIBO_MITO_PCA_GSEA) $(PCA_FLOW_STAMPS) \
	$(NO_RIBO_MITO_PCA_FLOW_STAMPS)

quinn_%.pca.eigenvectors.tsv quinn_%.pca.time: $(FILTERED_DIR)/quinn_%.tree.nex $(FILTERED_DIR)/quinn_%.moran.filtered.tsv
	/usr/bin/time -o quinn_$*.pca.time $(GEXPCA) \
		--tree $(FILTERED_DIR)/quinn_$*.tree.nex \
		--expr $(FILTERED_DIR)/quinn_$*.moran.filtered.tsv \
		--outprefix quinn_$* \
		--dim $(PCA_DIM) > quinn_$*.pca.term

quinn_%_no_ribo_mito.pca.eigenvectors.tsv quinn_%_no_ribo_mito.pca.time: \
$(FILTERED_DIR)/quinn_%.tree.nex $(FILTERED_DIR)/quinn_%_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o quinn_$*_no_ribo_mito.pca.time $(GEXPCA) \
		--tree $(FILTERED_DIR)/quinn_$*.tree.nex \
		--expr $(FILTERED_DIR)/quinn_$*_no_ribo_mito.moran.filtered.tsv \
		--outprefix quinn_$*_no_ribo_mito \
		--dim $(PCA_DIM) > quinn_$*_no_ribo_mito.pca.term

quinn_%.pca.gsea.tsv: quinn_%.pca.eigenvectors.tsv $(SRC)/run_factor_gsea.R $(CLUSTERPROFILER_SIF)
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/run_factor_gsea.R \
		/mnt/cwd/$< human /mnt/cwd/quinn_$*.pca.gsea \
		genes_by_components 10 \
		> quinn_$*.pca.gsea.term

quinn_%.pca.latent_flow.tsv: $(FILTERED_DIR)/quinn_%.moran.filtered.tsv quinn_%.pca.eigenvectors.tsv $(PREPROCESSED_DIR)/quinn_%.vine.mean.nwk $(SRC)/pca_scores_to_latent_flow.py $(HOTSPOT_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FILTERED_DIR):/mnt/gex_filter \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(HOTSPOT_SIF) \
		python /mnt/src/pca_scores_to_latent_flow.py \
		/mnt/gex_filter/quinn_$*.moran.filtered.tsv \
		/mnt/cwd/quinn_$*.pca.eigenvectors.tsv \
		/mnt/preprocessed_data/quinn_$*.vine.mean.nwk \
		/mnt/cwd/$@ \
		--max-components $(PCA_FLOW_DIMS)

quinn_%_no_ribo_mito.pca.latent_flow.tsv: $(FILTERED_DIR)/quinn_%_no_ribo_mito.moran.filtered.tsv quinn_%_no_ribo_mito.pca.eigenvectors.tsv $(PREPROCESSED_DIR)/quinn_%.vine.mean.nwk $(SRC)/pca_scores_to_latent_flow.py $(HOTSPOT_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FILTERED_DIR):/mnt/gex_filter \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(HOTSPOT_SIF) \
		python /mnt/src/pca_scores_to_latent_flow.py \
		/mnt/gex_filter/quinn_$*_no_ribo_mito.moran.filtered.tsv \
		/mnt/cwd/quinn_$*_no_ribo_mito.pca.eigenvectors.tsv \
		/mnt/preprocessed_data/quinn_$*.vine.mean.nwk \
		/mnt/cwd/$@ \
		--max-components $(PCA_FLOW_DIMS)

quinn_%.pca.latent_flow.plots.done: quinn_%.pca.latent_flow.tsv $(PREPROCESSED_DIR)/quinn_%_barcode_tissue.tsv $(SRC)/plot_latent_flow.py $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow.py \
		/mnt/cwd/$< \
		/mnt/cwd/quinn_$*.pca.latent_flow \
		--barcode-discrete-label-tsv /mnt/preprocessed_data/quinn_$*_barcode_tissue.tsv
	touch $@

quinn_%_no_ribo_mito.pca.latent_flow.plots.done: quinn_%_no_ribo_mito.pca.latent_flow.tsv $(PREPROCESSED_DIR)/quinn_%_barcode_tissue.tsv $(SRC)/plot_latent_flow.py $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow.py \
		/mnt/cwd/$< \
		/mnt/cwd/quinn_$*_no_ribo_mito.pca.latent_flow \
		--barcode-discrete-label-tsv /mnt/preprocessed_data/quinn_$*_barcode_tissue.tsv
	touch $@

clean:
	rm -f quinn_*.pca.*
