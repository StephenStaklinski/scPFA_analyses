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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/hotspot
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data

SAMPLES := quinn_1 quinn_2 quinn_3 quinn_4
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
CP1_MODULE_FLOW := quinn_1.hotspot.module_scores.latent_flow.tsv
CP1_MODULE_TREES := quinn_1.hotspot.module_trees.M3_M4_M5.600dpi.png

.PHONY: all archive clean cp1-module-trees
all: $(MODULES) $(GO_ENRICHMENT) $(SAMPLE_CORRELATIONS)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( *.hotspot.* module_sample_correlations.pdf ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No Hotspot outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

cp1-module-trees: $(CP1_MODULE_TREES)

$(CP1_MODULE_FLOW): \
		quinn_1.hotspot.module_scores.tsv \
		$(PREPROCESSED_DIR)/quinn_1.vine.mean.nwk \
		$(SRC)/hotspot_module_scores_to_latent_flow.py \
		$(HOTSPOT_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(PREPROCESSED_DIR):/mnt/preprocessed_data \
		$(HOTSPOT_SIF) \
		python /mnt/src/hotspot_module_scores_to_latent_flow.py \
		/mnt/cwd/quinn_1.hotspot.module_scores.tsv \
		/mnt/preprocessed_data/quinn_1.vine.mean.nwk \
		/mnt/cwd/$@

$(CP1_MODULE_TREES): $(CP1_MODULE_FLOW) $(SRC)/plot_latent_flow_trees.py $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow_trees.py \
		/mnt/cwd/$(CP1_MODULE_FLOW) \
		/mnt/cwd/quinn_1.hotspot.module_trees \
		--factor 3 --factor 4 --factor 5 \
		--title-prefix Module \
		--output /mnt/cwd/$@ \
		--dpi 600

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
		human \
		/mnt/cwd/$*.hotspot.modules.go \
		> $*.hotspot.modules.go.term 2>&1

clean:
	rm -f *.hotspot.*
