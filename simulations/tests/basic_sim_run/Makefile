export SHELL=/usr/bin/bash

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := ../data

all: sim.expr.tsv

# Main run on the example dataset
sim.term sim.time sim.expr.tsv: ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex
	/usr/bin/time -o sim.time \
		${GEX_LINEAGE_DIR}/bin/gexSim \
		--trees $< \
		--outprefix sim \
		--include-factorization > sim.term

# Remove the output files from the current run
clean:
	rm -f sim.*
