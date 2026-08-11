SHELL := /usr/bin/bash

MAIN_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/simulations/gene_phylo_filter_validation_sims
THREADS ?= 10

.PHONY: run-all collect-results collect-performance collect-time collect-treestats

run-all:
	for dir in $(MAIN_DIR)/*taxa; do \
		if [[ -d $$dir ]]; then \
			$(MAKE) -C $$dir -j $(THREADS) all; \
		fi; \
	done

run_treestats:
	for dir in $(MAIN_DIR)/*taxa; do \
		if [[ -d $$dir ]]; then \
			$(MAKE) -C $$dir -j $(THREADS) treestats; \
		fi; \
	done

collect-results: collect-performance collect-time collect-treestats

collect-performance:
	OUTFILE=$(MAIN_DIR)/eval.all.performance.txt; \
	rm -f $$OUTFILE; \
	performance_files=$$(find $(MAIN_DIR) -mindepth 2 -maxdepth 2 ! -path "$(MAIN_DIR)/archive*" -type f -name "eval.all.performance.txt"); \
	for file in $$performance_files; do \
		if [[ ! -f $$OUTFILE ]]; then \
			head -n 1 $$file > $$OUTFILE; \
		fi; \
		tail -n +2 $$file >> $$OUTFILE; \
	done; \
	python $(MAIN_DIR)/plot_all_performance.py $$OUTFILE $(MAIN_DIR)/eval.all.performance.pdf

collect-time:
	OUTFILE=$(MAIN_DIR)/eval.all.time.txt; \
	rm -f $$OUTFILE; \
	time_files=$$(find $(MAIN_DIR) -mindepth 2 -maxdepth 2 ! -path "$(MAIN_DIR)/archive*" -type f -name "eval.all.time.txt"); \
	for file in $$time_files; do \
		if [[ ! -f $$OUTFILE ]]; then \
			head -n 1 $$file > $$OUTFILE; \
		fi; \
		tail -n +2 $$file >> $$OUTFILE; \
	done; \
	python $(MAIN_DIR)/plot_all_time.py $$OUTFILE $(MAIN_DIR)/eval.all.time.pdf

collect-treestats: collect-performance
	OUTFILE=$(MAIN_DIR)/eval.all.tree_stats.tsv; \
	rm -f $$OUTFILE; \
	tree_stats_files=$$(find $(MAIN_DIR) -mindepth 2 -maxdepth 2 ! -path "$(MAIN_DIR)/archive*" -type f -name "eval.all.tree_stats.tsv"); \
	for file in $$tree_stats_files; do \
		if [[ ! -f $$OUTFILE ]]; then \
			head -n 1 $$file > $$OUTFILE; \
		fi; \
		tail -n +2 $$file >> $$OUTFILE; \
	done; \
	python $(MAIN_DIR)/plot_recall_pairwise_min.py \
		$(MAIN_DIR)/eval.all.performance.txt \
		$$OUTFILE \
		$(MAIN_DIR)/eval.all.tree_stats.pdf
