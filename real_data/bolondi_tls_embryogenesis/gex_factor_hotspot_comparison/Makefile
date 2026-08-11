export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif

CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_hotspot_comparison
FACTOR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_final
HOTSPOT_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/hotspot
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
HOTSPOT_SIF := $(CONTAINERS)/hotspot/hotspot.sif

PROGRAM_CORRELATIONS := gex_factor_hotspot_comparison.pdf
NO_RIBO_MITO_PROGRAM_CORRELATIONS := gex_factor_hotspot_comparison_no_ribo_mito.pdf
ACTIVE_NO_RIBO_MITO_PROGRAM_CORRELATIONS := gex_factor_hotspot_comparison_no_ribo_mito_active_factors.pdf

.PHONY: all residual-test archive clean
all: $(PROGRAM_CORRELATIONS) $(NO_RIBO_MITO_PROGRAM_CORRELATIONS) \
	$(ACTIVE_NO_RIBO_MITO_PROGRAM_CORRELATIONS)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( gex_factor_hotspot_comparison* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No gexFactor-Hotspot comparison outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

$(PROGRAM_CORRELATIONS): \
		$(HOTSPOT_DIR)/tls1.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls1.hotspot.autocorrelations.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.autocorrelations.tsv \
		$(FACTOR_DIR)/tls1.fit.L.tsv \
		$(FACTOR_DIR)/tls2.fit.L.tsv \
		$(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--sample TLS1 /mnt/hotspot/tls1.hotspot.modules.tsv /mnt/hotspot/tls1.hotspot.autocorrelations.tsv /mnt/factor/tls1.fit.L.tsv \
		--sample TLS2 /mnt/hotspot/tls2.hotspot.modules.tsv /mnt/hotspot/tls2.hotspot.autocorrelations.tsv /mnt/factor/tls2.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison

$(NO_RIBO_MITO_PROGRAM_CORRELATIONS): \
		$(HOTSPOT_DIR)/tls1.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls1.hotspot.autocorrelations.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.autocorrelations.tsv \
		$(FACTOR_DIR)/tls1_no_ribo_mito.fit.L.tsv \
		$(FACTOR_DIR)/tls2_no_ribo_mito.fit.L.tsv \
		$(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--sample TLS1_no_ribo_mito /mnt/hotspot/tls1.hotspot.modules.tsv /mnt/hotspot/tls1.hotspot.autocorrelations.tsv /mnt/factor/tls1_no_ribo_mito.fit.L.tsv \
		--sample TLS2_no_ribo_mito /mnt/hotspot/tls2.hotspot.modules.tsv /mnt/hotspot/tls2.hotspot.autocorrelations.tsv /mnt/factor/tls2_no_ribo_mito.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_no_ribo_mito

$(ACTIVE_NO_RIBO_MITO_PROGRAM_CORRELATIONS): \
		$(HOTSPOT_DIR)/tls1.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.modules.tsv \
		$(HOTSPOT_DIR)/tls1.hotspot.autocorrelations.tsv \
		$(HOTSPOT_DIR)/tls2.hotspot.autocorrelations.tsv \
		$(FACTOR_DIR)/tls1_no_ribo_mito.fit.L.tsv \
		$(FACTOR_DIR)/tls2_no_ribo_mito.fit.L.tsv \
		$(FACTOR_DIR)/tls1_no_ribo_mito.fit.summary.tsv \
		$(FACTOR_DIR)/tls2_no_ribo_mito.fit.summary.tsv \
		$(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--active-sample TLS1_no_ribo_mito /mnt/hotspot/tls1.hotspot.modules.tsv /mnt/hotspot/tls1.hotspot.autocorrelations.tsv /mnt/factor/tls1_no_ribo_mito.fit.L.tsv /mnt/factor/tls1_no_ribo_mito.fit.summary.tsv \
		--active-sample TLS2_no_ribo_mito /mnt/hotspot/tls2.hotspot.modules.tsv /mnt/hotspot/tls2.hotspot.autocorrelations.tsv /mnt/factor/tls2_no_ribo_mito.fit.L.tsv /mnt/factor/tls2_no_ribo_mito.fit.summary.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_no_ribo_mito_active_factors

clean:
	rm -f $(PROGRAM_CORRELATIONS) \
		gex_factor_hotspot_comparison.*.correlations.tsv \
		$(NO_RIBO_MITO_PROGRAM_CORRELATIONS) \
		gex_factor_hotspot_comparison_no_ribo_mito.*.correlations.tsv \
		$(ACTIVE_NO_RIBO_MITO_PROGRAM_CORRELATIONS) \
		gex_factor_hotspot_comparison_no_ribo_mito_active_factors.*.correlations.tsv \
		$(RESIDUAL_AUTOCORRELATIONS)
