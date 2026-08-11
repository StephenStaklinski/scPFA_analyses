export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
SCANPY_SIF := $(CONTAINERS)/scanpy/scanpy.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/pca_umap
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
CLONES := 1 2 3 4

PCA_PLOTS := $(CLONES:%=quinn_%_gex.pca.pdf)
UMAP_PLOTS := $(CLONES:%=quinn_%_gex.umap.pdf)

.PHONY: all archive clean
all: $(PCA_PLOTS) $(UMAP_PLOTS)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( quinn_*_gex.pca.pdf quinn_*_gex.pca.tsv quinn_*_gex.umap.pdf quinn_*_gex.umap.tsv ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No PCA or UMAP outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

quinn_%_gex.pca.pdf quinn_%_gex.umap.pdf quinn_%_gex.pca.tsv quinn_%_gex.umap.tsv: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/scanpy_dim_reduction.py /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

clean:
	rm -f quinn_*_gex.pca.pdf quinn_*_gex.pca.tsv quinn_*_gex.umap.pdf quinn_*_gex.umap.tsv
