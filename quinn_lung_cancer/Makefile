export SHELL=/usr/bin/bash

QUINNROOT := ../data

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:


BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
GEXLINEAGE := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexLineage

CURR_DIR := $(BENCHMARKS_DIR)/quinn_lung_cancer

DATA_DIR_CHAR := ./data/GSE161363/GSM4905334_character_matrices_5k
DATA_DIR_GEX := ./data/GSE161363/GSM4905335_barcodes_genes_matrix_5k


# Excluded in the original paper barcode quality control: 5 16 18 25 33 38 39 41 50 53 65 69 75 81 87 88 93
# Excluded because likely still too large due to VINE memory requirements: 1 2 3
QUINNCPS := 4 6 7 8 9 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 28 29 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100


CHARMATRICES := $(QUINNCPS:%=quinn_%_character_matrix.tsv)
GEXMATS := $(QUINNCPS:%=quinn_%_gex.tsv)
VINETREES := $(QUINNCPS:%=quinn.%.vine.nex)
GEXLINEAGERESULTS := $(QUINNCPS:%=quinn.%.fit.summary.tsv)


all: prep_data run_vine

prep_data: $(CHARMATRICES) $(GEXMATS)
run_vine: $(VINETREES)
gexlineage: $(GEXLINEAGERESULTS)

# Copy and reformat the provided character matrices
quinn_%_character_matrix.tsv:
	awk 'BEGIN{FS=OFS="\t"} { for(i=2;i<=NF;i++) if($$i=="-") $$i="-1"; print }' $(DATA_DIR_CHAR)/m5k_lg$*_character_matrix.alleleThresh.txt > $@

# Subset the gene expression data to the cells present in the character matrix
quinn_%_gex.tsv: quinn_%_character_matrix.tsv $(DATA_DIR_GEX)/GSM4905335_matrix.5k.mtx $(DATA_DIR_GEX)/GSM4905335_barcodes.5k.tsv $(DATA_DIR_GEX)/GSM4905335_genes.5k.tsv
		singularity exec --bind $(DATA_DIR_GEX):/mnt/gex --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/quinn_$*_character_matrix.tsv /mnt/gex/GSM4905335_matrix.5k.mtx \
			/mnt/gex/GSM4905335_barcodes.5k.tsv /mnt/gex/GSM4905335_genes.5k.tsv \
			/mnt/cwd/quinn_$*_gex.tsv > quinn_$*_gex.term

# Run VINE
quinn.%.vine.log quinn.%.vine.nwk quinn.%.vine.time quinn.%.vine.mean.nwk: quinn_%_character_matrix.tsv
	/usr/bin/time -o quinn.$*.vine.time $(VINE) -M 400 -c 100 -v 0 \
	 -i CRISPR quinn_$*_character_matrix.tsv --logf quinn.$*.vine.log --mean quinn.$*.vine.mean.nwk > quinn.$*.vine.nwk

# Convert nwk to nexus
quinn.%.vine.nex: quinn.%.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd/ --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py /mnt/cwd/quinn.$*.vine.nwk /mnt/cwd/quinn.$*.vine.nex

# Fit the gene expression model
quinn.%.fit.summary.tsv quinn.%.fit.F.tsv quinn.%.fit.L.tsv: quinn.%.vine.nex quinn_%_gex.tsv
	NTAXA=$$(wc -l < quinn_$*_gex.tsv); \
	if [ $$NTAXA -gt 1000 ]; then \
		FILTERTEST=moran; \
	else \
		FILTERTEST=lrt; \
	fi; \
	/usr/bin/time -o quinn.$*.fit.time ${GEXLINEAGE} \
	--trees quinn.$*.vine.nex \
	--expr quinn_$*_gex.tsv \
	--tree-total-time 1.0 \
	--outprefix quinn.$*.fit \
	--filter-test $$FILTERTEST \
	--seed 1 > quinn.$*.fit.term


