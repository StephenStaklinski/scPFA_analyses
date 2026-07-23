export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
DATA_DIR_CHAR := ../raw_data/GSE161363/GSM4905334_character_matrices_5k
DATA_DIR_GEX := ../raw_data/GSE161363/GSM4905335_barcodes_genes_matrix_5k
CLONES := 1 2 3 4

.PHONY: all
all: $(CLONES:%=quinn_%.vine.nex) $(CLONES:%=quinn_%_gex.tsv) $(CLONES:%=quinn_%_barcode_tissue.tsv)

quinn_%_character_matrix.tsv:
	awk 'BEGIN{FS=OFS="\t"} { for(i=2;i<=NF;i++) if($$i=="-") $$i="-1"; print }' $(DATA_DIR_CHAR)/m5k_lg$*_character_matrix.alleleThresh.txt > $@

quinn_%_barcode_tissue.tsv: quinn_%_character_matrix.tsv
	awk 'BEGIN{FS=OFS="\t"; print "barcode", "tissue"} FNR > 1 { tissue=$$1; sub(/\..*$$/, "", tissue); print $$1, tissue }' $< > $@

quinn_%_gex.tsv: quinn_%_character_matrix.tsv $(DATA_DIR_GEX)/GSM4905335_matrix.5k.mtx $(DATA_DIR_GEX)/GSM4905335_barcodes.5k.tsv $(DATA_DIR_GEX)/GSM4905335_genes.5k.tsv
	singularity exec --bind $(abspath $(DATA_DIR_GEX)):/mnt/gex --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(CASSIOPEIA_SIF) python /mnt/src/subset_mtx_to_char_matrix_cells.py /mnt/cwd/quinn_$*_character_matrix.tsv /mnt/gex/GSM4905335_matrix.5k.mtx /mnt/gex/GSM4905335_barcodes.5k.tsv /mnt/gex/GSM4905335_genes.5k.tsv /mnt/cwd/$@ > quinn_$*_gex.term

quinn_%.vine.nwk: quinn_%_character_matrix.tsv
	/usr/bin/time -o quinn_$*.vine.time $(VINE) -M 400 -c 100 -v 0 -i CRISPR $< --logf quinn_$*.vine.log --mean quinn_$*.vine.mean.nwk > $@

quinn_%.vine.nex: quinn_%.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) python /mnt/src/nwk2nex.py /mnt/cwd/$< /mnt/cwd/$@
