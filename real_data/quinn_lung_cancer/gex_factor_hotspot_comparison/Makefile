export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif

CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_hotspot_comparison
FACTOR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_final
HOTSPOT_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/hotspot
SAMPLES := quinn_1 quinn_2 quinn_3 quinn_4

OUTPUT := gex_factor_hotspot_comparison.pdf
NO_RIBO_MITO_OUTPUT := gex_factor_hotspot_comparison_no_ribo_mito.pdf
ACTIVE_NO_RIBO_MITO_OUTPUT := gex_factor_hotspot_comparison_no_ribo_mito_active_factors.pdf
CP1_CP2_ACTIVE_NO_RIBO_MITO_OUTPUT := gex_factor_hotspot_comparison_cp1_cp2_no_ribo_mito_active_factors.pdf
CP1_CP2_CP4_ACTIVE_NO_RIBO_MITO_OUTPUT := gex_factor_hotspot_comparison_cp1_cp2_cp4_no_ribo_mito_active_factors.pdf
HOTSPOT_INPUTS := $(foreach sample,$(SAMPLES),\
	$(HOTSPOT_DIR)/$(sample).hotspot.modules.tsv \
	$(HOTSPOT_DIR)/$(sample).hotspot.autocorrelations.tsv)
FACTOR_INPUTS := $(foreach sample,$(SAMPLES),\
	$(FACTOR_DIR)/$(sample).fit.L.tsv)
NO_RIBO_MITO_FACTOR_INPUTS := $(foreach sample,$(SAMPLES),\
	$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.L.tsv)
NO_RIBO_MITO_SUMMARY_INPUTS := $(foreach sample,$(SAMPLES),\
	$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.summary.tsv)

.PHONY: all archive clean

all: $(OUTPUT) $(NO_RIBO_MITO_OUTPUT) $(ACTIVE_NO_RIBO_MITO_OUTPUT) \
	$(CP1_CP2_ACTIVE_NO_RIBO_MITO_OUTPUT) \
	$(CP1_CP2_CP4_ACTIVE_NO_RIBO_MITO_OUTPUT)

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

$(OUTPUT): $(HOTSPOT_INPUTS) $(FACTOR_INPUTS) $(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--sample CP1 /mnt/hotspot/quinn_1.hotspot.modules.tsv /mnt/hotspot/quinn_1.hotspot.autocorrelations.tsv /mnt/factor/quinn_1.fit.L.tsv \
		--sample CP2 /mnt/hotspot/quinn_2.hotspot.modules.tsv /mnt/hotspot/quinn_2.hotspot.autocorrelations.tsv /mnt/factor/quinn_2.fit.L.tsv \
		--sample CP3 /mnt/hotspot/quinn_3.hotspot.modules.tsv /mnt/hotspot/quinn_3.hotspot.autocorrelations.tsv /mnt/factor/quinn_3.fit.L.tsv \
		--sample CP4 /mnt/hotspot/quinn_4.hotspot.modules.tsv /mnt/hotspot/quinn_4.hotspot.autocorrelations.tsv /mnt/factor/quinn_4.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison

$(NO_RIBO_MITO_OUTPUT): $(HOTSPOT_INPUTS) $(NO_RIBO_MITO_FACTOR_INPUTS) $(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--sample CP1_no_ribo_mito /mnt/hotspot/quinn_1.hotspot.modules.tsv /mnt/hotspot/quinn_1.hotspot.autocorrelations.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.L.tsv \
		--sample CP2_no_ribo_mito /mnt/hotspot/quinn_2.hotspot.modules.tsv /mnt/hotspot/quinn_2.hotspot.autocorrelations.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.L.tsv \
		--sample CP3_no_ribo_mito /mnt/hotspot/quinn_3.hotspot.modules.tsv /mnt/hotspot/quinn_3.hotspot.autocorrelations.tsv /mnt/factor/quinn_3_no_ribo_mito.fit.L.tsv \
		--sample CP4_no_ribo_mito /mnt/hotspot/quinn_4.hotspot.modules.tsv /mnt/hotspot/quinn_4.hotspot.autocorrelations.tsv /mnt/factor/quinn_4_no_ribo_mito.fit.L.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_no_ribo_mito

$(ACTIVE_NO_RIBO_MITO_OUTPUT): $(HOTSPOT_INPUTS) $(NO_RIBO_MITO_FACTOR_INPUTS) \
		$(NO_RIBO_MITO_SUMMARY_INPUTS) $(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--active-sample CP1_no_ribo_mito /mnt/hotspot/quinn_1.hotspot.modules.tsv /mnt/hotspot/quinn_1.hotspot.autocorrelations.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.summary.tsv \
		--active-sample CP2_no_ribo_mito /mnt/hotspot/quinn_2.hotspot.modules.tsv /mnt/hotspot/quinn_2.hotspot.autocorrelations.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.summary.tsv \
		--active-sample CP3_no_ribo_mito /mnt/hotspot/quinn_3.hotspot.modules.tsv /mnt/hotspot/quinn_3.hotspot.autocorrelations.tsv /mnt/factor/quinn_3_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_3_no_ribo_mito.fit.summary.tsv \
		--active-sample CP4_no_ribo_mito /mnt/hotspot/quinn_4.hotspot.modules.tsv /mnt/hotspot/quinn_4.hotspot.autocorrelations.tsv /mnt/factor/quinn_4_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_4_no_ribo_mito.fit.summary.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_no_ribo_mito_active_factors

$(CP1_CP2_ACTIVE_NO_RIBO_MITO_OUTPUT): \
		$(foreach sample,quinn_1 quinn_2,$(HOTSPOT_DIR)/$(sample).hotspot.modules.tsv \
			$(HOTSPOT_DIR)/$(sample).hotspot.autocorrelations.tsv \
			$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.L.tsv \
			$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.summary.tsv) \
		$(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--active-sample CP1_no_ribo_mito /mnt/hotspot/quinn_1.hotspot.modules.tsv /mnt/hotspot/quinn_1.hotspot.autocorrelations.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.summary.tsv \
		--active-sample CP2_no_ribo_mito /mnt/hotspot/quinn_2.hotspot.modules.tsv /mnt/hotspot/quinn_2.hotspot.autocorrelations.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.summary.tsv \
		--columns 2 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_cp1_cp2_no_ribo_mito_active_factors

$(CP1_CP2_CP4_ACTIVE_NO_RIBO_MITO_OUTPUT): \
		$(foreach sample,quinn_1 quinn_2 quinn_4,$(HOTSPOT_DIR)/$(sample).hotspot.modules.tsv \
			$(HOTSPOT_DIR)/$(sample).hotspot.autocorrelations.tsv \
			$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.L.tsv \
			$(FACTOR_DIR)/$(sample)_no_ribo_mito.fit.summary.tsv) \
		$(SRC)/plot_factor_hotspot_comparison.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		--bind $(HOTSPOT_DIR):/mnt/hotspot \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_hotspot_comparison.py \
		--active-sample CP1_no_ribo_mito /mnt/hotspot/quinn_1.hotspot.modules.tsv /mnt/hotspot/quinn_1.hotspot.autocorrelations.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.summary.tsv \
		--active-sample CP2_no_ribo_mito /mnt/hotspot/quinn_2.hotspot.modules.tsv /mnt/hotspot/quinn_2.hotspot.autocorrelations.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.summary.tsv \
		--active-sample CP4_no_ribo_mito /mnt/hotspot/quinn_4.hotspot.modules.tsv /mnt/hotspot/quinn_4.hotspot.autocorrelations.tsv /mnt/factor/quinn_4_no_ribo_mito.fit.L.tsv /mnt/factor/quinn_4_no_ribo_mito.fit.summary.tsv \
		--columns 3 \
		--out-prefix /mnt/cwd/gex_factor_hotspot_comparison_cp1_cp2_cp4_no_ribo_mito_active_factors

clean:
	rm -f $(OUTPUT) gex_factor_hotspot_comparison.*.correlations.tsv \
		$(NO_RIBO_MITO_OUTPUT) \
		gex_factor_hotspot_comparison_no_ribo_mito.*.correlations.tsv \
		$(ACTIVE_NO_RIBO_MITO_OUTPUT) \
		gex_factor_hotspot_comparison_no_ribo_mito_active_factors.*.correlations.tsv \
		$(CP1_CP2_ACTIVE_NO_RIBO_MITO_OUTPUT) \
		gex_factor_hotspot_comparison_cp1_cp2_no_ribo_mito_active_factors.*.correlations.tsv \
		$(CP1_CP2_CP4_ACTIVE_NO_RIBO_MITO_OUTPUT) \
		gex_factor_hotspot_comparison_cp1_cp2_cp4_no_ribo_mito_active_factors.*.correlations.tsv \
		gex_factor_hotspot_comparison_cp1_cp4_no_ribo_mito_active_factors.pdf \
		gex_factor_hotspot_comparison_cp1_cp4_no_ribo_mito_active_factors.*.correlations.tsv \
		gex_factor_hotspot_comparison_cp2_cp4_no_ribo_mito_active_factors.pdf \
		gex_factor_hotspot_comparison_cp2_cp4_no_ribo_mito_active_factors.*.correlations.tsv
