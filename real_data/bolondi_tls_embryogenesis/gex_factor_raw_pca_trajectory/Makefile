export SHELL := /usr/bin/bash

.DELETE_ON_ERROR:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
PYPLOTTING_SIF := $(BENCHMARKS_DIR)/containers/pyplotting/pyplotting.sif
PCA_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/pca_umap
FIT_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_final

SAMPLES := tls1 tls2
SCRIPT := plot_factor_distance_raw_pca.py
OUTPUTS := $(SAMPLES:%=%.factor_latent_distance.raw_pca.pdf)

.PHONY: all clean

all: $(OUTPUTS)

%.factor_latent_distance.raw_pca.pdf: \
		$(PCA_DIR)/%_gex.pca.tsv \
		$(FIT_DIR)/%_no_ribo_mito.fit.latent_flow.tsv \
		$(FIT_DIR)/%_no_ribo_mito.fit.summary.tsv \
		$(SCRIPT) $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURDIR):/mnt/cwd \
		--bind $(PCA_DIR):/mnt/pca:ro \
		--bind $(FIT_DIR):/mnt/fit:ro \
		$(PYPLOTTING_SIF) \
		python /mnt/cwd/$(SCRIPT) \
			/mnt/pca/$*_gex.pca.tsv \
			/mnt/fit/$*_no_ribo_mito.fit.latent_flow.tsv \
			/mnt/fit/$*_no_ribo_mito.fit.summary.tsv \
			/mnt/cwd/$@

clean:
	rm -f $(OUTPUTS)
