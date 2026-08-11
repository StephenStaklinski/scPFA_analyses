export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
DATA_DIR := ../raw_data/GSE220949
TLS1_PATH := $(DATA_DIR)/GSM6831757_TracerTLS1
TLS2_PATH := $(DATA_DIR)/GSM6831759_TracerTLS2
SAMPLES := tls1 tls2

.PHONY: all archive
all: $(SAMPLES:%=%.vine.nex) $(SAMPLES:%=%_gex.tsv)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( \
		tls*_character_matrix.tsv \
		tls*_mutation_dict.tsv \
		tls*_gex.tsv \
		tls*_character_matrix.term \
		tls*_gex.term \
		tls*.vine.* \
	); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No preprocessing or VINE outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

tls1_character_matrix.tsv tls1_mutation_dict.tsv: $(TLS1_PATH)/allele_table.txt
	singularity exec --bind $(abspath $<):/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(CASSIOPEIA_SIF) python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls1_character_matrix.tsv /mnt/cwd/tls1_mutation_dict.tsv > tls1_character_matrix.term

tls2_character_matrix.tsv tls2_mutation_dict.tsv: $(TLS2_PATH)/allele_table.txt
	singularity exec --bind $(abspath $<):/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(CASSIOPEIA_SIF) python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls2_character_matrix.tsv /mnt/cwd/tls2_mutation_dict.tsv > tls2_character_matrix.term

tls1_gex.tsv: $(TLS1_PATH)/matrix.mtx $(TLS1_PATH)/barcodes.tsv $(TLS1_PATH)/features.tsv tls1_character_matrix.tsv
	singularity exec --bind $(abspath $(TLS1_PATH)):/mnt/input --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(CASSIOPEIA_SIF) python /mnt/src/subset_mtx_to_char_matrix_cells.py /mnt/cwd/tls1_character_matrix.tsv /mnt/input/matrix.mtx /mnt/input/barcodes.tsv /mnt/input/features.tsv /mnt/cwd/$@ > tls1_gex.term

tls2_gex.tsv: $(TLS2_PATH)/matrix.mtx $(TLS2_PATH)/barcodes.tsv $(TLS2_PATH)/features.tsv tls2_character_matrix.tsv
	singularity exec --bind $(abspath $(TLS2_PATH)):/mnt/input --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(CASSIOPEIA_SIF) python /mnt/src/subset_mtx_to_char_matrix_cells.py /mnt/cwd/tls2_character_matrix.tsv /mnt/input/matrix.mtx /mnt/input/barcodes.tsv /mnt/input/features.tsv /mnt/cwd/$@ > tls2_gex.term

%.vine.nwk: %_character_matrix.tsv
	/usr/bin/time -o $*.vine.time $(VINE) -M 400 -c 100 -v 0 -i CRISPR $< --logf $*.vine.log --mean $*.vine.mean.nwk > $@

%.vine.nex: %.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) python /mnt/src/nwk2nex.py /mnt/cwd/$< /mnt/cwd/$@
