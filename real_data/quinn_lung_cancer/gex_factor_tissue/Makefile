export SHELL := /usr/bin/bash

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_filter
FACTOR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_final
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_tissue

TREE_OUTPUT := cp1_cp2.tissue_colored_trees.pdf
ASSOCIATION_OUTPUT := cp1_cp2.active_factor_tissue_association.pdf
ASSOCIATION_TABLE := cp1_cp2.active_factor_tissue_association.tsv
CLASSIFIER_OUTPUT := cp1_cp2.active_factor_tissue_classifier.pdf
CLASSIFIER_METRICS := cp1_cp2.active_factor_tissue_classifier.metrics.tsv
CLASSIFIER_PREDICTIONS := cp1_cp2.active_factor_tissue_classifier.predictions.tsv
TREE_INPUTS := $(FILTERED_DIR)/quinn_1.tree.nex $(FILTERED_DIR)/quinn_2.tree.nex
ASSOCIATION_INPUTS := \
	$(FACTOR_DIR)/quinn_1_no_ribo_mito.fit.F.tsv \
	$(FACTOR_DIR)/quinn_1_no_ribo_mito.fit.summary.tsv \
	$(FACTOR_DIR)/quinn_2_no_ribo_mito.fit.F.tsv \
	$(FACTOR_DIR)/quinn_2_no_ribo_mito.fit.summary.tsv

.PHONY: all archive clean

all: $(TREE_OUTPUT) $(ASSOCIATION_OUTPUT) $(CLASSIFIER_OUTPUT)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( cp1_cp2.tissue_colored_trees.* cp1_cp2.active_factor_tissue_* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No tissue-analysis outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

$(TREE_OUTPUT): $(TREE_INPUTS) $(SRC)/plot_tissue_colored_trees.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FILTERED_DIR):/mnt/trees \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_tissue_colored_trees.py \
		--tree CP1 /mnt/trees/quinn_1.tree.nex \
		--tree CP2 /mnt/trees/quinn_2.tree.nex \
		--output /mnt/cwd/$@

$(ASSOCIATION_OUTPUT): $(ASSOCIATION_INPUTS) $(SRC)/plot_factor_tissue_association.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_tissue_association.py \
		--sample CP1 /mnt/factor/quinn_1_no_ribo_mito.fit.F.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.summary.tsv \
		--sample CP2 /mnt/factor/quinn_2_no_ribo_mito.fit.F.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.summary.tsv \
		--output /mnt/cwd/$(ASSOCIATION_OUTPUT) \
		--table-output /mnt/cwd/$(ASSOCIATION_TABLE)

$(ASSOCIATION_TABLE): $(ASSOCIATION_OUTPUT)
	@test -e $@

$(CLASSIFIER_OUTPUT): $(ASSOCIATION_INPUTS) $(SRC)/classify_tissue_from_factors.py
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		--bind $(FACTOR_DIR):/mnt/factor \
		$(PYPLOTTING_SIF) \
		python /mnt/src/classify_tissue_from_factors.py \
		--sample CP1 /mnt/factor/quinn_1_no_ribo_mito.fit.F.tsv /mnt/factor/quinn_1_no_ribo_mito.fit.summary.tsv \
		--sample CP2 /mnt/factor/quinn_2_no_ribo_mito.fit.F.tsv /mnt/factor/quinn_2_no_ribo_mito.fit.summary.tsv \
		--output /mnt/cwd/$(CLASSIFIER_OUTPUT) \
		--metrics-output /mnt/cwd/$(CLASSIFIER_METRICS) \
		--predictions-output /mnt/cwd/$(CLASSIFIER_PREDICTIONS)

$(CLASSIFIER_METRICS) $(CLASSIFIER_PREDICTIONS): $(CLASSIFIER_OUTPUT)
	@test -e $@

clean:
	rm -f $(TREE_OUTPUT) $(ASSOCIATION_OUTPUT) $(ASSOCIATION_TABLE) \
		$(CLASSIFIER_OUTPUT) $(CLASSIFIER_METRICS) $(CLASSIFIER_PREDICTIONS)
