export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
GEXLINEAGE := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexLineage

CURR_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis

DATA_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis/data/GSE220949
TLS1_PATH := $(DATA_DIR)/GSM6831757_TracerTLS1
TLS2_PATH := $(DATA_DIR)/GSM6831759_TracerTLS2
TLSMULTI_PATH := $(DATA_DIR)/GSM6831760_TracerTLSMulti

CHARMATRICES := tls1_character_matrix.tsv tls2_character_matrix.tsv
GEXMATS := tls1_gex.tsv tls2_gex.tsv
VINETREES := tls1.vine.nex tls2.vine.nex
GEXLINEAGERESULTS := tls1.fit.summary.tsv tls2.fit.summary.tsv

all: prep_data run_vine

prep_data: $(CHARMATRICES) $(GEXMATS)
run_vine: $(VINETREES)
gexlineage: $(GEXLINEAGERESULTS)

# Obtain the TLS1 character matrix
tls1_character_matrix.tsv tls1_mutation_dict.tsv: $(TLS1_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls1_character_matrix.tsv /mnt/cwd/tls1_mutation_dict.tsv > tls1_character_matrix.term

# Obtain the TLS2 character matrix
tls2_character_matrix.tsv tls2_mutation_dict.tsv: $(TLS2_PATH)/allele_table.txt
		singularity exec --bind $^:/mnt/allele_table.txt --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/cassiopeia_allele_table_to_character_matrix.py /mnt/allele_table.txt /mnt/cwd/tls2_character_matrix.tsv /mnt/cwd/tls2_mutation_dict.tsv > tls2_character_matrix.term

# TODO: Split the TLSMulti allele table into separate tables for each of the additional 21 samples, obtain the character matrices, and fill in the other steps below

# Subset the TLS1 gene expression data to the cells present in the character matrix
tls1_gex.tsv: $(TLS1_PATH)/matrix.mtx $(TLS1_PATH)/barcodes.tsv $(TLS1_PATH)/features.tsv tls1_character_matrix.tsv
		singularity exec --bind $(TLS1_PATH):/mnt/tls1 --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/tls1_character_matrix.tsv /mnt/tls1/matrix.mtx \
			/mnt/tls1/barcodes.tsv /mnt/tls1/features.tsv /mnt/cwd/tls1_gex.tsv > tls1_gex.term

# Subset the TLS2 gene expression data to the cells present in the character matrix
tls2_gex.tsv: $(TLS2_PATH)/matrix.mtx $(TLS2_PATH)/barcodes.tsv $(TLS2_PATH)/features.tsv tls2_character_matrix.tsv
		singularity exec --bind $(TLS2_PATH):/mnt/tls2 --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/tls2_character_matrix.tsv /mnt/tls2/matrix.mtx \
			/mnt/tls2/barcodes.tsv /mnt/tls2/features.tsv /mnt/cwd/tls2_gex.tsv > tls2_gex.term

# Run VINE on the character matrices to infer lineage trees
%.vine.nwk: %_character_matrix.tsv
	/usr/bin/time -o $*.vine.time $(VINE) -M 400 -c 100 -v 0 \
	-i CRISPR $*_character_matrix.tsv --logf $*.vine.log --mean $*.vine.mean.nwk > $*.vine.nwk

# Convert nwk to nexus
%.vine.nex: %.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd/ --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py /mnt/cwd/$*.vine.nwk /mnt/cwd/$*.vine.nex

# Fit the gene expression model
%.fit.summary.tsv %.fit.F.tsv %.fit.L.tsv: %.vine.nex %_gex.tsv
	/usr/bin/time -o $*.fit.time ${GEXLINEAGE} \
	--trees $*.vine.nex \
	--expr $*_gex.tsv \
	--tree-total-time 1.0 \
	--outprefix $*.fit \
	--filter-test moran \
	--seed 1 > $*.fit.term

