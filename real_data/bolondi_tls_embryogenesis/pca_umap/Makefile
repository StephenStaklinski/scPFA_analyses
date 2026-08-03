export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
SCANPY_SIF := $(CONTAINERS)/scanpy/scanpy.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/pca_umap
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
SAMPLES := tls1 tls2

PCA_PLOTS := $(SAMPLES:%=%_gex.pca.pdf)
UMAP_PLOTS := $(SAMPLES:%=%_gex.umap.pdf)

.PHONY: all clean
all: $(PCA_PLOTS) $(UMAP_PLOTS)

%_gex.pca.pdf %_gex.umap.pdf %_gex.pca.tsv %_gex.umap.tsv: $(PREPROCESSED_DIR)/%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/scanpy_dim_reduction.py /mnt/preprocessed_data/$*_gex.tsv /mnt/cwd/$*_gex

clean:
	rm -f *_gex.pca.pdf *_gex.pca.tsv *_gex.umap.pdf *_gex.umap.tsv
