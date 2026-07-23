export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
SCANPY_SIF := $(CONTAINERS)/scanpy/scanpy.sif
MONOCLE_SIF := $(CONTAINERS)/monocle/monocle.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/trajectory_velocity
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
SAMPLES := tls1 tls2

.PHONY: all palantir monocle
all: palantir monocle
palantir: $(SAMPLES:%=%_gex.palantir.pseudotime.tsv)
monocle: $(SAMPLES:%=%_gex.monocle.pseudotime.pdf)

%_gex.palantir.pseudotime.tsv: $(PREPROCESSED_DIR)/%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/palantir_trajectory.py /mnt/preprocessed_data/$*_gex.tsv /mnt/cwd/$*_gex

%_gex.monocle.pseudotime.pdf: $(PREPROCESSED_DIR)/%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(MONOCLE_SIF) Rscript /mnt/src/monocle_trajectory.R /mnt/preprocessed_data/$*_gex.tsv /mnt/cwd/$*_gex
