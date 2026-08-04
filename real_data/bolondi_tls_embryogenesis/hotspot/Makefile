export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
HOTSPOT_SIF := $(CONTAINERS)/hotspot/hotspot.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/hotspot
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data

SAMPLES := tls1 tls2
JOBS ?= 4
MIN_CELLS ?= 10
N_NEIGHBORS ?= 30
GENE_FDR ?= 0.05
MODULE_FDR ?= 0.05
MIN_MODULE_GENES ?= 50
MAX_GENES ?= 0

MODULES := $(SAMPLES:%=%.hotspot.modules.tsv)
GO_ENRICHMENT := $(SAMPLES:%=%.hotspot.modules.go.tsv)
SAMPLE_CORRELATIONS := module_sample_correlations.pdf

.PHONY: all clean
all: $(MODULES) $(GO_ENRICHMENT) $(SAMPLE_CORRELATIONS)

$(SAMPLE_CORRELATIONS): \
		$(MODULES) \
		$(SAMPLES:%=%.hotspot.autocorrelations.tsv) \
		$(SRC)/plot_hotspot_module_sample_correlations.py \
		$(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_hotspot_module_sample_correlations.py \
		$(foreach sample,$(SAMPLES),--sample $(sample) /mnt/cwd/$(sample).hotspot.modules.tsv /mnt/cwd/$(sample).hotspot.autocorrelations.tsv) \
		--pair-only \
		--output /mnt/cwd/$@

%.hotspot.modules.tsv: \
		$(PREPROCESSED_DIR)/%_gex.tsv \
		$(PREPROCESSED_DIR)/%.vine.mean.nwk \
		$(HOTSPOT_SIF)
	/usr/bin/time -o $*.hotspot.time \
		singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(HOTSPOT_SIF) \
		python /mnt/src/run_hotspot_tree.py \
		--counts /mnt/preprocessed_data/$*_gex.tsv \
		--tree /mnt/preprocessed_data/$*.vine.mean.nwk \
		--out-prefix /mnt/cwd/$*.hotspot \
		--jobs $(JOBS) \
		--min-cells $(MIN_CELLS) \
		--n-neighbors $(N_NEIGHBORS) \
		--gene-fdr $(GENE_FDR) \
		--module-fdr $(MODULE_FDR) \
		--min-module-genes $(MIN_MODULE_GENES) \
		--max-genes $(MAX_GENES) \
		> $*.hotspot.term 2>&1

%.hotspot.modules.go.tsv: %.hotspot.modules.tsv %.hotspot.autocorrelations.tsv
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/run_hotspot_module_go.R \
		/mnt/cwd/$*.hotspot.modules.tsv \
		/mnt/cwd/$*.hotspot.autocorrelations.tsv \
		mouse \
		/mnt/cwd/$*.hotspot.modules.go \
		> $*.hotspot.modules.go.term 2>&1

clean:
	rm -f *.hotspot.*
