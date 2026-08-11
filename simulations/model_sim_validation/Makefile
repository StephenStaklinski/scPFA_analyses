SHELL := /usr/bin/bash

MAIN_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/simulations/model_sim_validation
THREADS ?= 10

.PHONY: all diff times archive-all run-all rerun-eval restore-archive

run-all:
	for dir in $(MAIN_DIR)/*taxa*genes; do \
		if [[ -d $$dir ]]; then \
			$(MAKE) -C $$dir -j $(THREADS) all; \
		fi; \
	done

rerun-eval:
	for dir in $(MAIN_DIR)/*taxa*genes; do \
		if [[ -d $$dir ]]; then \
			rm -f $$dir/sim.*.eval.summary.tsv $$dir/sim.*.eval.term $$dir/diff.txt; \
			$(MAKE) -C $$dir -j $(THREADS) eval; \
		fi; \
	done

restore-archive:
	if [[ -z "$(ARCHIVE)" ]]; then \
		echo "Usage: make restore-archive ARCHIVE=/path/to/archive"; \
		exit 1; \
	fi; \
	archive="$(ARCHIVE)"; \
	if [[ ! -d $$archive ]]; then \
		echo "Archive directory not found: $$archive"; \
		exit 1; \
	fi; \
	shopt -s nullglob; \
	top_files=( \
		"$$archive"/diff.txt \
		"$$archive"/diff.pdf \
		"$$archive"/eval.all.time.txt \
	); \
	if (( $${#top_files[@]} > 0 )); then \
		mv "$${top_files[@]}" $(MAIN_DIR)/; \
	fi; \
	for archived_dir in "$$archive"/*taxa*genes; do \
		if [[ -d $$archived_dir ]]; then \
			dest="$(MAIN_DIR)/$$(basename "$$archived_dir")"; \
			mkdir -p $$dest; \
			files=( "$$archived_dir"/sim* "$$archived_dir"/eval.all* "$$archived_dir"/diff.txt ); \
			if (( $${#files[@]} > 0 )); then \
				mv "$${files[@]}" "$$dest"/; \
			fi; \
		fi; \
	done; \
	shopt -u nullglob

collect-results: diff

diff: times
	OUTFILE=$(MAIN_DIR)/diff.txt; \
	rm -f $$OUTFILE; \
	diff_files=$$(find $(MAIN_DIR) -mindepth 2 -maxdepth 2 ! -path "$(MAIN_DIR)/archive*" -type f -name "diff.txt"); \
	for file in $$diff_files; do \
		if [[ ! -f $$OUTFILE ]]; then \
			head -n 1 $$file > $$OUTFILE; \
		fi; \
		tail -n +2 $$file >> $$OUTFILE; \
	done; \
	python $(MAIN_DIR)/plot_diff_time.py \
		$$OUTFILE \
		$(MAIN_DIR)/eval.all.time.txt \
		$(MAIN_DIR)/diff.pdf

times:
	OUTFILE=$(MAIN_DIR)/eval.all.time.txt; \
	rm -f $$OUTFILE; \
	time_files=$$(find $(MAIN_DIR) -mindepth 2 -maxdepth 2 ! -path "$(MAIN_DIR)/archive*" -type f -name "eval.all.time.txt"); \
	for file in $$time_files; do \
		if [[ ! -f $$OUTFILE ]]; then \
			head -n 1 $$file > $$OUTFILE; \
		fi; \
		tail -n +2 $$file >> $$OUTFILE; \
	done

archive-all:
	archive_dir=$(MAIN_DIR)/archive_all_$$(date +%Y-%m-%d_%H-%M); \
	mkdir -p $$archive_dir; \
	for file in \
		$(MAIN_DIR)/diff.txt \
		$(MAIN_DIR)/diff.pdf \
		$(MAIN_DIR)/eval.all.time.txt; do \
		if [[ -e $$file ]]; then \
			mv $$file $$archive_dir/; \
		fi; \
	done; \
	for dir in $(MAIN_DIR)/*taxa*genes; do \
		if [[ -d $$dir ]]; then \
			subdir=$$archive_dir/$$(basename $$dir); \
			mkdir -p $$subdir; \
			shopt -s nullglob; \
			files=( "$$dir"/sim* "$$dir"/eval.all* "$$dir"/diff.txt ); \
			if (( $${#files[@]} > 0 )); then \
				mv "$${files[@]}" $$subdir/; \
			fi; \
			shopt -u nullglob; \
		fi; \
	done

clean-all:
	rm -f $(MAIN_DIR)/diff.txt $(MAIN_DIR)/diff.pdf $(MAIN_DIR)/eval.all.time.txt
	for dir in $(MAIN_DIR)/*taxa*genes; do \
		if [[ -d $$dir ]]; then \
			make -C $$dir clean; \
		fi; \
	done
