export SHELL := /bin/bash

BENCHMARK_DIR := $(abspath ../../..)
SRC_DIR := $(BENCHMARK_DIR)/src
PCA_DIR := $(abspath ../gex_pca)
FACTOR_DIR := $(abspath ../gex_factor_final)
CURR_DIR := $(abspath .)
PYPLOTTING_SIF := $(BENCHMARK_DIR)/containers/pyplotting/pyplotting.sif

OUTPUT := gex_pca_gex_factor_comparison.pdf
NO_RIBO_MITO_OUTPUT := gex_pca_gex_factor_comparison_no_ribo_mito.pdf
PLOT_SCRIPT := $(SRC_DIR)/plot_pca_factor_comparison.py
SAMPLES := tls1 tls2
PCA_INPUTS := $(foreach sample,$(SAMPLES),\
	$(PCA_DIR)/$(sample).pca.eigenvectors.tsv)
FACTOR_INPUTS := $(foreach sample,$(SAMPLES),\
	$(FACTOR_DIR)/$(sample).fit.L.tsv)
NO_RIBO_MITO_PCA_INPUTS := $(foreach sample,$(SAMPLES),\
	$(PCA_DIR)/$(sample)_no_ribo_mito.pca.eigenvectors.tsv)
NO_RIBO_MITO_FACTOR_INPUTS := $(foreach sample,$(SAMPLES),\
	$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.L.tsv)

.PHONY: all archive clean

all: $(OUTPUT) $(NO_RIBO_MITO_OUTPUT)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( gex_pca_gex_factor_comparison* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No PCA-gexFactor comparison outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

$(OUTPUT): $(PCA_INPUTS) $(FACTOR_INPUTS) $(PLOT_SCRIPT) $(PYPLOTTING_SIF)
	/usr/bin/time -o gex_pca_gex_factor_comparison.time \
		singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(PCA_DIR):/mnt/pca \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(SRC_DIR):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_pca_factor_comparison.py \
		--sample TLS1 /mnt/pca/tls1.pca.eigenvectors.tsv /mnt/factor/tls1.fit.L.tsv \
		--sample TLS2 /mnt/pca/tls2.pca.eigenvectors.tsv /mnt/factor/tls2.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_pca_gex_factor_comparison \
		> gex_pca_gex_factor_comparison.term 2>&1

$(NO_RIBO_MITO_OUTPUT): $(NO_RIBO_MITO_PCA_INPUTS) $(NO_RIBO_MITO_FACTOR_INPUTS) $(PLOT_SCRIPT) $(PYPLOTTING_SIF)
	/usr/bin/time -o gex_pca_gex_factor_comparison_no_ribo_mito.time \
		singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(PCA_DIR):/mnt/pca \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(SRC_DIR):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_pca_factor_comparison.py \
		--sample TLS1_no_ribo_mito /mnt/pca/tls1_no_ribo_mito.pca.eigenvectors.tsv /mnt/factor/tls1_no_ribo_mito.fit.L.tsv \
		--sample TLS2_no_ribo_mito /mnt/pca/tls2_no_ribo_mito.pca.eigenvectors.tsv /mnt/factor/tls2_no_ribo_mito.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_pca_gex_factor_comparison_no_ribo_mito \
		> gex_pca_gex_factor_comparison_no_ribo_mito.term 2>&1

clean:
	rm -f gex_pca_gex_factor_comparison.pdf \
		gex_pca_gex_factor_comparison.*.correlations.tsv \
		gex_pca_gex_factor_comparison.time \
		gex_pca_gex_factor_comparison.term \
		gex_pca_gex_factor_comparison_no_ribo_mito.pdf \
		gex_pca_gex_factor_comparison_no_ribo_mito.*.correlations.tsv \
		gex_pca_gex_factor_comparison_no_ribo_mito.time \
		gex_pca_gex_factor_comparison_no_ribo_mito.term
