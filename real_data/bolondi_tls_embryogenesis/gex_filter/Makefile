export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXFILTER := $(GEX_LINEAGE_DIR)/bin/gexFilter
BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_filter
PREPROCESSED_DIR := ../preprocessed_data

SAMPLES := tls1 tls2
TREE_INDEX ?= 1
FILTER_TEST ?= moran
SEED ?= 1

TREES := $(SAMPLES:%=%.tree.nex)
FILTERED_GEX := $(SAMPLES:%=%.$(FILTER_TEST).filtered.tsv)
NO_RIBO_MITO_FILTERED_GEX := $(SAMPLES:%=%_no_ribo_mito.$(FILTER_TEST).filtered.tsv)
NO_RIBO_MITO_MORAN_RESULTS := $(SAMPLES:%=%_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv)
RIBO_MITO_HERITABILITY_PLOT := ribo_mito_heritability.pdf
GO_PLOTS := \
	no_ribo_mito.heritable.go.pdf \
	no_ribo_mito.not_heritable.go.pdf
GO_SAMPLE_RESULTS := $(SAMPLES:%=%_no_ribo_mito.mouse.go.tsv)
GO_RESULTS := no_ribo_mito.go.tsv
OVERLAP_PLOTS := no_ribo_mito.heritable_overlap.pdf

.PHONY: all filter filter_no_ribo_mito go overlap ribo_mito_heritability archive clean
all: filter go overlap ribo_mito_heritability
filter: $(FILTERED_GEX) filter_no_ribo_mito
filter_no_ribo_mito: $(NO_RIBO_MITO_FILTERED_GEX)
go: $(GO_PLOTS)
overlap: $(OVERLAP_PLOTS)
ribo_mito_heritability: $(RIBO_MITO_HERITABILITY_PLOT)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( \
		*.tree.nex *.filtered.tsv *.correlation.*.tsv *.term *.time \
		*.go.pdf *.go.tsv *.heritable_overlap.pdf *.heritable_overlap.tsv \
		ribo_mito_heritability.pdf \
	); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No gexFilter outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

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

%_no_ribo_mito.$(FILTER_TEST).filtered.tsv: \
%.tree.nex $(PREPROCESSED_DIR)/%_gex.tsv
	prefix=$(@:.filtered.tsv=); \
	/usr/bin/time -o $$prefix.time $(GEXFILTER) \
		--trees $< \
		--expr $(PREPROCESSED_DIR)/$*_gex.tsv \
		--outprefix $$prefix \
		--filter-test $(FILTER_TEST) \
		--remove-ribo-mito-genes \
		--seed $(SEED) > $$prefix.term

%.$(FILTER_TEST).correlation.moran.tsv: %.$(FILTER_TEST).filtered.tsv
	@test -e $@ || { echo "ERROR: expected Moran results missing: $@" >&2; exit 1; }

%_no_ribo_mito.mouse.go.tsv: \
%_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv $(CLUSTERPROFILER_SIF)
	sample_label=$*; \
	sample_label=$${sample_label^^}; \
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/calculate_filter_go_enrichment.R \
		both /mnt/cwd/$*_no_ribo_mito.mouse.go \
		mouse \
		$$sample_label \
		/mnt/cwd/$*_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		> $*_no_ribo_mito.mouse.go.calculate.term 2>&1

no_ribo_mito.go.tsv: $(GO_SAMPLE_RESULTS)
	awk 'FNR == 1 && NR != 1 { next } { print }' $^ > $@

no_ribo_mito.heritable.go.pdf: $(GO_RESULTS) $(CLUSTERPROFILER_SIF)
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/plot_filter_go_enrichment.R \
		/mnt/cwd/no_ribo_mito.go.tsv \
		/mnt/cwd/no_ribo_mito.heritable.go.pdf \
		passing \
		TLS1 TLS2 \
		> no_ribo_mito.heritable.go.plot.term 2>&1

no_ribo_mito.not_heritable.go.pdf: $(GO_RESULTS) $(CLUSTERPROFILER_SIF)
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/plot_filter_go_enrichment.R \
		/mnt/cwd/no_ribo_mito.go.tsv \
		/mnt/cwd/no_ribo_mito.not_heritable.go.pdf \
		not_passing \
		TLS1 TLS2 \
		> no_ribo_mito.not_heritable.go.plot.term 2>&1

no_ribo_mito.heritable_overlap.pdf: $(NO_RIBO_MITO_MORAN_RESULTS) $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_filter_gene_set_overlap.py \
		/mnt/cwd/tls1_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/tls2_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/no_ribo_mito.heritable_overlap \
		--sample-a TLS1 \
		--sample-b TLS2

ribo_mito_heritability.pdf: \
$(SAMPLES:%=%.$(FILTER_TEST).correlation.moran.tsv) \
$(NO_RIBO_MITO_MORAN_RESULTS) \
$(SRC)/plot_ribo_mito_heritability.py $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_ribo_mito_heritability.py \
		/mnt/cwd/$@ \
		--sample TLS1 \
			/mnt/cwd/tls1.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/tls1_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		--sample TLS2 \
			/mnt/cwd/tls2.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/tls2_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv

clean:
	rm -f *.tree.nex *.filtered.tsv *.correlation.*.tsv *.term *.time
	rm -f *.go.pdf *.go.tsv
	rm -f *.heritable_overlap.pdf *.heritable_overlap.tsv
	rm -f *.ribo_mito_heritability.pdf
