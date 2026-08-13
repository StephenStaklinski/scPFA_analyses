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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_filter
PREPROCESSED_DIR := ../preprocessed_data

CLONES := 1 2 3 4
FILTER_TEST ?= moran
SEED ?= 1

FILTERED_GEX := $(CLONES:%=quinn_%.$(FILTER_TEST).filtered.tsv)
NO_RIBO_MITO_FILTERED_GEX := $(CLONES:%=quinn_%_no_ribo_mito.$(FILTER_TEST).filtered.tsv)
MORAN_RESULTS := $(CLONES:%=quinn_%.$(FILTER_TEST).correlation.moran.tsv)
NO_RIBO_MITO_MORAN_RESULTS := $(CLONES:%=quinn_%_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv)
RIBO_MITO_HERITABILITY_PLOT := ribo_mito_heritability.pdf
GO_RESULTS := no_ribo_mito.go.tsv
GO_SAMPLE_RESULTS := $(CLONES:%=quinn_%_no_ribo_mito.human.go.tsv)
GO_PLOTS := \
	no_ribo_mito.heritable.go.pdf \
	no_ribo_mito.not_heritable.go.pdf
OVERLAP_PLOTS := no_ribo_mito.heritable_overlap.pdf
OVERLAP_SCRIPT := $(SRC)/plot_filter_gene_set_overlap_four_samples.py

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
		quinn_*.filtered.tsv quinn_*.correlation.*.tsv \
		quinn_*.term quinn_*.time *.go.pdf *.go.tsv \
		*.heritable_overlap.pdf *.heritable_overlap.tsv ribo_mito_heritability.pdf \
	); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No gexFilter outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

quinn_%.$(FILTER_TEST).filtered.tsv: \
$(PREPROCESSED_DIR)/quinn_%.vine.mean.nex $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	prefix=$(@:.filtered.tsv=); \
	/usr/bin/time -o $$prefix.time $(GEXFILTER) \
		--trees $< \
		--expr $(PREPROCESSED_DIR)/quinn_$*_gex.tsv \
		--outprefix $$prefix \
		--filter-test $(FILTER_TEST) \
		--seed $(SEED) > $$prefix.term

quinn_%_no_ribo_mito.$(FILTER_TEST).filtered.tsv: \
$(PREPROCESSED_DIR)/quinn_%.vine.mean.nex $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	prefix=$(@:.filtered.tsv=); \
	/usr/bin/time -o $$prefix.time $(GEXFILTER) \
		--trees $< \
		--expr $(PREPROCESSED_DIR)/quinn_$*_gex.tsv \
		--outprefix $$prefix \
		--filter-test $(FILTER_TEST) \
		--remove-ribo-mito-genes \
		--seed $(SEED) > $$prefix.term

quinn_%.$(FILTER_TEST).correlation.moran.tsv: \
quinn_%.$(FILTER_TEST).filtered.tsv
	@test -e $@ || { echo "ERROR: expected Moran results missing: $@" >&2; exit 1; }

quinn_%_no_ribo_mito.human.go.tsv: \
quinn_%_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
$(CLUSTERPROFILER_SIF)
	singularity run \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(CLUSTERPROFILER_SIF) \
		/mnt/src/calculate_filter_go_enrichment.R \
		both /mnt/cwd/quinn_$*_no_ribo_mito.human.go \
		human \
		CP$* \
		/mnt/cwd/quinn_$*_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		> quinn_$*_no_ribo_mito.human.go.calculate.term 2>&1

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
		CP1 CP2 CP3 CP4 \
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
		CP1 CP2 CP3 CP4 \
		> no_ribo_mito.not_heritable.go.plot.term 2>&1

no_ribo_mito.heritable_overlap.pdf: \
$(NO_RIBO_MITO_FILTERED_GEX) $(PYPLOTTING_SIF)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_filter_gene_set_overlap_four_samples.py \
		/mnt/cwd/quinn_1_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/quinn_2_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/quinn_3_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/quinn_4_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		/mnt/cwd/no_ribo_mito.heritable_overlap \
		--sample-names CP1 CP2 CP3 CP4

ribo_mito_heritability.pdf: \
$(MORAN_RESULTS) $(NO_RIBO_MITO_MORAN_RESULTS)
	singularity exec \
		--bind $(CURR_DIR):/mnt/cwd \
		--bind $(SRC):/mnt/src \
		$(PYPLOTTING_SIF) \
		python /mnt/src/plot_ribo_mito_heritability.py \
		/mnt/cwd/$@ \
		--sample CP1 \
			/mnt/cwd/quinn_1.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/quinn_1_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		--sample CP2 \
			/mnt/cwd/quinn_2.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/quinn_2_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		--sample CP3 \
			/mnt/cwd/quinn_3.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/quinn_3_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv \
		--sample CP4 \
			/mnt/cwd/quinn_4.$(FILTER_TEST).correlation.moran.tsv \
			/mnt/cwd/quinn_4_no_ribo_mito.$(FILTER_TEST).correlation.moran.tsv

clean:
	rm -f quinn_*.filtered.tsv \
		quinn_*.correlation.*.tsv quinn_*.term quinn_*.time
	rm -f *.go.pdf *.go.tsv
	rm -f *.heritable_overlap.pdf *.heritable_overlap.tsv
	rm -f ribo_mito_heritability.pdf *go.plot.term
