export SHELL := /usr/bin/bash

.DELETE_ON_ERROR:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
PYPLOTTING_SIF := $(BENCHMARKS_DIR)/containers/pyplotting/pyplotting.sif
PCA_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/pca_umap
FIT_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_final

CLONES := 1 2 3 4
SCRIPT := plot_factor_distance_raw_pca.py
OUTPUTS := $(CLONES:%=quinn_%.factor_latent_distance.raw_pca.pdf)

.PHONY: all clean

all: $(OUTPUTS)

quinn_%.factor_latent_distance.raw_pca.pdf: \
		$(PCA_DIR)/quinn_%_gex.pca.tsv \
		$(FIT_DIR)/quinn_%_no_ribo_mito.fit.latent_flow.tsv \
		$(FIT_DIR)/quinn_%_no_ribo_mito.fit.summary.tsv \
		$(SCRIPT) $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURDIR):/mnt/cwd \
		--bind $(PCA_DIR):/mnt/pca:ro \
		--bind $(FIT_DIR):/mnt/fit:ro \
		$(PYPLOTTING_SIF) \
		python /mnt/cwd/$(SCRIPT) \
			/mnt/pca/quinn_$*_gex.pca.tsv \
			/mnt/fit/quinn_$*_no_ribo_mito.fit.latent_flow.tsv \
			/mnt/fit/quinn_$*_no_ribo_mito.fit.summary.tsv \
			/mnt/cwd/$@

clean:
	rm -f $(OUTPUTS)
