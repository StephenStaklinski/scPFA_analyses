export SHELL := /usr/bin/bash

.DELETE_ON_ERROR:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
PYPLOTTING_SIF := $(BENCHMARKS_DIR)/containers/pyplotting/pyplotting.sif
TLS_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_final
CP_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_final

SCRIPT := plot_cross_context_factor_comparison.py
TREE_SCRIPT := plot_active_factor_ancestral_trees.py
HOTSPOT_SCRIPT := plot_cross_context_hotspot_module_comparison.py
OUT_PREFIX := bolondi_quinn_active_factors
TREE_OUTPUT := bolondi_quinn_active_factor_ancestral_trees.png
HOTSPOT_PREFIX := bolondi_quinn_hotspot_modules

INPUTS := \
	$(TLS_DIR)/tls1_no_ribo_mito.fit.L.gsea.tsv \
	$(TLS_DIR)/tls1_no_ribo_mito.fit.summary.tsv \
	$(TLS_DIR)/tls2_no_ribo_mito.fit.L.gsea.tsv \
	$(TLS_DIR)/tls2_no_ribo_mito.fit.summary.tsv \
	$(CP_DIR)/quinn_1_no_ribo_mito.fit.L.gsea.tsv \
	$(CP_DIR)/quinn_1_no_ribo_mito.fit.summary.tsv \
	$(CP_DIR)/quinn_2_no_ribo_mito.fit.L.gsea.tsv \
	$(CP_DIR)/quinn_2_no_ribo_mito.fit.summary.tsv

OUTPUTS := \
	$(OUT_PREFIX).cross_context_factor_similarity.pdf \
	$(OUT_PREFIX).cross_context_factor_similarity.tsv

HOTSPOT_OUTPUTS := \
	$(HOTSPOT_PREFIX).cross_context_hotspot_module_similarity.pdf \
	$(HOTSPOT_PREFIX).cross_context_hotspot_module_similarity.tsv

.PHONY: all clean

all: $(OUTPUTS) $(TREE_OUTPUT) $(HOTSPOT_OUTPUTS)

$(OUTPUTS): $(SCRIPT) $(INPUTS) $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURDIR):/mnt/cwd \
		--bind $(TLS_DIR):/mnt/tls:ro \
		--bind $(CP_DIR):/mnt/cp:ro \
		$(PYPLOTTING_SIF) \
		python /mnt/cwd/$(SCRIPT) \
			--sample TLS1 development /mnt/tls/tls1_no_ribo_mito.fit.L.gsea.tsv /mnt/tls/tls1_no_ribo_mito.fit.summary.tsv \
			--sample TLS2 development /mnt/tls/tls2_no_ribo_mito.fit.L.gsea.tsv /mnt/tls/tls2_no_ribo_mito.fit.summary.tsv \
			--sample CP1 cancer /mnt/cp/quinn_1_no_ribo_mito.fit.L.gsea.tsv /mnt/cp/quinn_1_no_ribo_mito.fit.summary.tsv \
			--sample CP2 cancer /mnt/cp/quinn_2_no_ribo_mito.fit.L.gsea.tsv /mnt/cp/quinn_2_no_ribo_mito.fit.summary.tsv \
			--out-prefix /mnt/cwd/$(OUT_PREFIX)

$(HOTSPOT_OUTPUTS): $(HOTSPOT_SCRIPT) \
		$(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/hotspot/tls1.hotspot.modules.go.tsv \
		$(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/hotspot/tls2.hotspot.modules.go.tsv \
		$(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/hotspot/quinn_1.hotspot.modules.go.tsv \
		$(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/hotspot/quinn_2.hotspot.modules.go.tsv \
		$(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURDIR):/mnt/cwd \
		--bind $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/hotspot:/mnt/tls_hotspot:ro \
		--bind $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/hotspot:/mnt/cp_hotspot:ro \
		$(PYPLOTTING_SIF) \
		python /mnt/cwd/$(HOTSPOT_SCRIPT) \
			--sample TLS1 development /mnt/tls_hotspot/tls1.hotspot.modules.go.tsv \
			--sample TLS2 development /mnt/tls_hotspot/tls2.hotspot.modules.go.tsv \
			--sample CP1 cancer /mnt/cp_hotspot/quinn_1.hotspot.modules.go.tsv \
			--sample CP2 cancer /mnt/cp_hotspot/quinn_2.hotspot.modules.go.tsv \
			--out-prefix /mnt/cwd/$(HOTSPOT_PREFIX)

$(TREE_OUTPUT): $(TREE_SCRIPT) \
		$(TLS_DIR)/tls1_no_ribo_mito.fit.latent_flow.tsv $(TLS_DIR)/tls1_no_ribo_mito.fit.summary.tsv \
		$(TLS_DIR)/tls2_no_ribo_mito.fit.latent_flow.tsv $(TLS_DIR)/tls2_no_ribo_mito.fit.summary.tsv \
		$(CP_DIR)/quinn_1_no_ribo_mito.fit.latent_flow.tsv $(CP_DIR)/quinn_1_no_ribo_mito.fit.summary.tsv \
		$(CP_DIR)/quinn_2_no_ribo_mito.fit.latent_flow.tsv $(CP_DIR)/quinn_2_no_ribo_mito.fit.summary.tsv \
		$(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURDIR):/mnt/cwd \
		--bind $(TLS_DIR):/mnt/tls:ro \
		--bind $(CP_DIR):/mnt/cp:ro \
		$(PYPLOTTING_SIF) \
		python /mnt/cwd/$(TREE_SCRIPT) \
			--clone TLS1 /mnt/tls/tls1_no_ribo_mito.fit.latent_flow.tsv /mnt/tls/tls1_no_ribo_mito.fit.summary.tsv \
			--clone TLS2 /mnt/tls/tls2_no_ribo_mito.fit.latent_flow.tsv /mnt/tls/tls2_no_ribo_mito.fit.summary.tsv \
			--clone CP1 /mnt/cp/quinn_1_no_ribo_mito.fit.latent_flow.tsv /mnt/cp/quinn_1_no_ribo_mito.fit.summary.tsv \
			--clone CP2 /mnt/cp/quinn_2_no_ribo_mito.fit.latent_flow.tsv /mnt/cp/quinn_2_no_ribo_mito.fit.summary.tsv \
			--output /mnt/cwd/$@

clean:
	rm -f $(OUTPUTS) $(TREE_OUTPUT) $(HOTSPOT_OUTPUTS)
