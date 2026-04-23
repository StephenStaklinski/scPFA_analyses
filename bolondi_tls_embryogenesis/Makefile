export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif

CURR_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis

DATA_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis/data/GSE220949
TLS1_PATH := $(DATA_DIR)/GSM6831757_TracerTLS1
TLS2_PATH := $(DATA_DIR)/GSM6831759_TracerTLS2
TLSMULTI_PATH := $(DATA_DIR)/GSM6831760_TracerTLSMulti


# Obtain the TLS1 character matrix
tls1_character_matrix.tsv tls1_mutation_dict.tsv: $(TLS1_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls1_character_matrix.tsv /mnt/cwd/tls1_mutation_dict.tsv > tls1_character_matrix.term

# Obtain the TLS2 character matrix
tls2_character_matrix.tsv tls2_mutation_dict.tsv: $(TLS2_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls2_character_matrix.tsv /mnt/cwd/tls2_mutation_dict.tsv > tls2_character_matrix.term

