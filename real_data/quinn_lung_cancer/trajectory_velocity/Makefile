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

.PHONY: all palantir monocle
all: palantir monocle
palantir: $(CLONES:%=quinn_%_gex.palantir.pseudotime.tsv)
monocle: $(CLONES:%=quinn_%_gex.monocle.pseudotime.pdf)

quinn_%_gex.palantir.pseudotime.tsv: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/palantir_trajectory.py /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

quinn_%_gex.monocle.pseudotime.pdf: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(MONOCLE_SIF) Rscript /mnt/src/monocle_trajectory.R /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex
