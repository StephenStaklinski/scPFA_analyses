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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_pca
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_filter

CLONES := 1 2 3 4
PCA_DIM := 50
PCA_EIGENVECTORS := $(CLONES:%=quinn_%.pca.eigenvectors.tsv)
PCA_TIMES := $(CLONES:%=quinn_%.pca.time)
PCA_GSEA := $(CLONES:%=quinn_%.pca.gsea.tsv)
NO_RIBO_MITO_PCA_EIGENVECTORS := $(CLONES:%=quinn_%_no_ribo_mito.pca.eigenvectors.tsv)
NO_RIBO_MITO_PCA_TIMES := $(CLONES:%=quinn_%_no_ribo_mito.pca.time)
NO_RIBO_MITO_PCA_GSEA := $(CLONES:%=quinn_%_no_ribo_mito.pca.gsea.tsv)

.PHONY: all clean
all: $(PCA_EIGENVECTORS) $(PCA_TIMES) $(PCA_GSEA) \
	$(NO_RIBO_MITO_PCA_EIGENVECTORS) $(NO_RIBO_MITO_PCA_TIMES) \
	$(NO_RIBO_MITO_PCA_GSEA)

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

clean:
	rm -f quinn_*.pca.*
