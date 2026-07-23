export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXFILTER := $(GEX_LINEAGE_DIR)/bin/gexFilter
PREPROCESSED_DIR := ../preprocessed_data

SAMPLES := tls1 tls2
TREE_INDEX ?= 1
FILTER_TEST ?= moran
SEED ?= 1

TREES := $(SAMPLES:%=%.tree.nex)
FILTERED_GEX := $(SAMPLES:%=%.$(FILTER_TEST).filtered.tsv)

.PHONY: all filter clean
all: filter
filter: $(FILTERED_GEX)

%.tree.nex: $(PREPROCESSED_DIR)/%.vine.nex
	awk -v tree_index="$(TREE_INDEX)" '\
		BEGIN { tree_count = 0; found = 0; } \
		/^[[:space:]]*[Tt][Rr][Ee][Ee][[:space:]]/ { \
			tree_count++; \
			if (tree_count == tree_index) { print; found = 1; } \
			next; \
		} \
		{ print; } \
		END { \
			if (!found) { \
				printf "TREE_INDEX %d was not found in %s\n", tree_index, FILENAME > "/dev/stderr"; \
				exit 1; \
			} \
		}' $< > $@

%.$(FILTER_TEST).filtered.tsv: %.tree.nex $(PREPROCESSED_DIR)/%_gex.tsv
	prefix=$(@:.filtered.tsv=); \
	/usr/bin/time -o $$prefix.time $(GEXFILTER) \
		--trees $< \
		--expr $(PREPROCESSED_DIR)/$*_gex.tsv \
		--outprefix $$prefix \
		--filter-test $(FILTER_TEST) \
		--seed $(SEED) > $$prefix.term

clean:
	rm -f *.tree.nex *.filtered.tsv *.correlation.*.tsv *.term *.time
