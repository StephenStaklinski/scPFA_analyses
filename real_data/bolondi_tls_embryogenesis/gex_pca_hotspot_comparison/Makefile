export SHELL := /bin/bash

BENCHMARK_DIR := $(abspath ../../..)
SRC_DIR := $(BENCHMARK_DIR)/src
HOTSPOT_DIR := $(abspath ../hotspot)
PCA_DIR := $(abspath ../gex_pca)
CURR_DIR := $(abspath .)
PYPLOTTING_SIF := $(BENCHMARK_DIR)/containers/pyplotting/pyplotting.sif

OUTPUT := pca_hotspot_comparison.pdf
NO_RIBO_MITO_OUTPUT := pca_hotspot_comparison_no_ribo_mito.pdf
PLOT_SCRIPT := $(SRC_DIR)/plot_pca_hotspot_comparison.py
HOTSPOT_INPUTS := \
	$(HOTSPOT_DIR)/tls1.hotspot.modules.tsv \
	$(HOTSPOT_DIR)/tls1.hotspot.autocorrelations.tsv \
	$(HOTSPOT_DIR)/tls2.hotspot.modules.tsv \
	$(HOTSPOT_DIR)/tls2.hotspot.autocorrelations.tsv
PCA_INPUTS := \
	$(PCA_DIR)/tls1.pca.eigenvectors.tsv \
	$(PCA_DIR)/tls2.pca.eigenvectors.tsv
NO_RIBO_MITO_PCA_INPUTS := \
	$(PCA_DIR)/tls1_no_ribo_mito.pca.eigenvectors.tsv \
	$(PCA_DIR)/tls2_no_ribo_mito.pca.eigenvectors.tsv

.PHONY: all clean

all: $(OUTPUT) $(NO_RIBO_MITO_OUTPUT)

$(OUTPUT): $(HOTSPOT_INPUTS) $(PCA_INPUTS) $(PLOT_SCRIPT) $(PYPLOTTING_SIF)
	/usr/bin/time -o pca_hotspot_comparison.time \
		singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		--bind $(PCA_DIR):/mnt/pca \
		--bind $(SRC_DIR):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_pca_hotspot_comparison.py \
		--sample TLS1 /mnt/hotspot/tls1.hotspot.modules.tsv /mnt/hotspot/tls1.hotspot.autocorrelations.tsv /mnt/pca/tls1.pca.eigenvectors.tsv \
		--sample TLS2 /mnt/hotspot/tls2.hotspot.modules.tsv /mnt/hotspot/tls2.hotspot.autocorrelations.tsv /mnt/pca/tls2.pca.eigenvectors.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/pca_hotspot_comparison \
		> pca_hotspot_comparison.term 2>&1

$(NO_RIBO_MITO_OUTPUT): $(HOTSPOT_INPUTS) $(NO_RIBO_MITO_PCA_INPUTS) $(PLOT_SCRIPT) $(PYPLOTTING_SIF)
	/usr/bin/time -o pca_hotspot_comparison_no_ribo_mito.time \
		singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		--bind $(PCA_DIR):/mnt/pca \
		--bind $(SRC_DIR):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_pca_hotspot_comparison.py \
		--sample TLS1_no_ribo_mito /mnt/hotspot/tls1.hotspot.modules.tsv /mnt/hotspot/tls1.hotspot.autocorrelations.tsv /mnt/pca/tls1_no_ribo_mito.pca.eigenvectors.tsv \
		--sample TLS2_no_ribo_mito /mnt/hotspot/tls2.hotspot.modules.tsv /mnt/hotspot/tls2.hotspot.autocorrelations.tsv /mnt/pca/tls2_no_ribo_mito.pca.eigenvectors.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/pca_hotspot_comparison_no_ribo_mito \
		> pca_hotspot_comparison_no_ribo_mito.term 2>&1

clean:
	rm -f pca_hotspot_comparison.pdf \
		pca_hotspot_comparison.*.correlations.tsv \
		pca_hotspot_comparison.time \
		pca_hotspot_comparison.term \
		pca_hotspot_comparison_no_ribo_mito.pdf \
		pca_hotspot_comparison_no_ribo_mito.*.correlations.tsv \
		pca_hotspot_comparison_no_ribo_mito.time \
		pca_hotspot_comparison_no_ribo_mito.term
