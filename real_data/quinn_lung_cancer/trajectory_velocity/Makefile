export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
SCANPY_SIF := $(CONTAINERS)/scanpy/scanpy.sif
MONOCLE_SIF := $(CONTAINERS)/monocle/monocle.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/trajectory_velocity
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
CLONES := 1 2 3 4

.PHONY: all palantir monocle archive
all: palantir monocle
palantir: $(CLONES:%=quinn_%_gex.palantir.pseudotime.tsv)
monocle: $(CLONES:%=quinn_%_gex.monocle.pseudotime.pdf)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( quinn_*_gex.palantir* quinn_*_gex.monocle* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No trajectory or velocity outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

quinn_%_gex.palantir.pseudotime.tsv: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/palantir_trajectory.py /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

quinn_%_gex.monocle.pseudotime.pdf: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(MONOCLE_SIF) Rscript /mnt/src/monocle_trajectory.R /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex
