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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_pca
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_filter

SAMPLES := tls1 tls2
PCA_DIM := 50
PCA_EIGENVECTORS := $(SAMPLES:%=%.pca.eigenvectors.tsv)
PCA_TIMES := $(SAMPLES:%=%.pca.time)
PCA_GSEA := $(SAMPLES:%=%.pca.gsea.tsv)
NO_RIBO_MITO_SAMPLES := $(SAMPLES:%=%_no_ribo_mito)
NO_RIBO_MITO_PCA_EIGENVECTORS := $(NO_RIBO_MITO_SAMPLES:%=%.pca.eigenvectors.tsv)
NO_RIBO_MITO_PCA_TIMES := $(NO_RIBO_MITO_SAMPLES:%=%.pca.time)
NO_RIBO_MITO_PCA_GSEA := $(NO_RIBO_MITO_SAMPLES:%=%.pca.gsea.tsv)

.PHONY: all archive clean
all: $(PCA_EIGENVECTORS) $(PCA_TIMES) $(PCA_GSEA) \
	$(NO_RIBO_MITO_PCA_EIGENVECTORS) $(NO_RIBO_MITO_PCA_TIMES) \
	$(NO_RIBO_MITO_PCA_GSEA)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( *.pca.* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No PCA outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

%.pca.eigenvectors.tsv %.pca.time: $(PREPROCESSED_DIR)/%.vine.mean.nex $(FILTERED_DIR)/%.moran.filtered.tsv
	/usr/bin/time -o $*.pca.time $(GEXPCA) \
		--tree $(PREPROCESSED_DIR)/$*.vine.mean.nex \
		--expr $(FILTERED_DIR)/$*.moran.filtered.tsv \
		--outprefix $* \
		--dim $(PCA_DIM) > $*.pca.term

%_no_ribo_mito.pca.eigenvectors.tsv %_no_ribo_mito.pca.time: \
$(PREPROCESSED_DIR)/%.vine.mean.nex $(FILTERED_DIR)/%_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o $*_no_ribo_mito.pca.time $(GEXPCA) \
		--tree $(PREPROCESSED_DIR)/$*.vine.mean.nex \
		--expr $(FILTERED_DIR)/$*_no_ribo_mito.moran.filtered.tsv \
		--outprefix $*_no_ribo_mito \
		--dim $(PCA_DIM) > $*_no_ribo_mito.pca.term

%.pca.gsea.tsv: %.pca.eigenvectors.tsv $(SRC)/run_factor_gsea.R $(CLUSTERPROFILER_SIF)
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/run_factor_gsea.R \
		/mnt/cwd/$< mouse /mnt/cwd/$*.pca.gsea \
		genes_by_components 10 \
		> $*.pca.gsea.term

clean:
	rm -f *.pca.*
