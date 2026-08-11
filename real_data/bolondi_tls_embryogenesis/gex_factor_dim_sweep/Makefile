export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
GEXFACTOR := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexFactor
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_dim_sweep
PREPROCESSED_DIR := ../preprocessed_data
FILTERED_DIR := ../gex_filter

SAMPLES := tls1 tls2
NO_RIBO_MITO_SAMPLES := $(SAMPLES:%=%_no_ribo_mito)
DIM_SUMMARIES := $(SAMPLES:%=%.dim_summary.tsv) $(NO_RIBO_MITO_SAMPLES:%=%.dim_summary.tsv)
DIM_VALUES ?= 1 2 3 4 5 10 15
NO_RIBO_MITO_DIM_VALUES ?= $(DIM_VALUES) 20 25 30
DIM_FIT_SUMMARIES = $(foreach dim,$(2),$(1).K$(dim).fit.summary.tsv)
DIM_FIT_GSEAS = $(foreach dim,$(2),$(1).K$(dim).fit.L.gsea.tsv)
DIM_FIT_TIMES = $(foreach dim,$(2),$(1).K$(dim).fit.time)

.PHONY: all archive clean-gexlineage
all: tls.dim_summary.pdf tls_no_ribo_mito.dim_summary.pdf

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( *.K*.fit.* *.dim_summary.* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No dimension-sweep outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

tls1.K%.fit.summary.tsv: $(FILTERED_DIR)/tls1.tree.nex $(FILTERED_DIR)/tls1.moran.filtered.tsv
	/usr/bin/time -o tls1.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls1.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls1.K$*.fit.term

tls2.K%.fit.summary.tsv: $(FILTERED_DIR)/tls2.tree.nex $(FILTERED_DIR)/tls2.moran.filtered.tsv
	/usr/bin/time -o tls2.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls2.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls2.K$*.fit.term

tls1_no_ribo_mito.K%.fit.summary.tsv: $(FILTERED_DIR)/tls1.tree.nex $(FILTERED_DIR)/tls1_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o tls1_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls1_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls1_no_ribo_mito.K$*.fit.term

tls2_no_ribo_mito.K%.fit.summary.tsv: $(FILTERED_DIR)/tls2.tree.nex $(FILTERED_DIR)/tls2_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o tls2_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls2_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls2_no_ribo_mito.K$*.fit.term

tls1.K%.fit.L.tsv tls1.K%.fit.log tls1.K%.fit.time: tls1.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }
tls2.K%.fit.L.tsv tls2.K%.fit.log tls2.K%.fit.time: tls2.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }
tls1_no_ribo_mito.K%.fit.L.tsv tls1_no_ribo_mito.K%.fit.log tls1_no_ribo_mito.K%.fit.time: tls1_no_ribo_mito.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }
tls2_no_ribo_mito.K%.fit.L.tsv tls2_no_ribo_mito.K%.fit.log tls2_no_ribo_mito.K%.fit.time: tls2_no_ribo_mito.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

%.fit.L.gsea.tsv: %.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< mouse /mnt/cwd/$*.fit.L.gsea > $*.fit.gsea.term

define DIM_SUMMARY_RULE
$(1).dim_summary.tsv: $(call DIM_FIT_SUMMARIES,$(1),$(2)) $(call DIM_FIT_GSEAS,$(1),$(2)) \
		$(call DIM_FIT_TIMES,$(1),$(2)) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(2),$(1).K$(dim).fit.log); do test -e $$$$log || { echo "ERROR: expected fit log missing: $$$$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_dim_log_likelihood.py --summary-only \
		/mnt/cwd/$(1).dim_summary.pdf \
		$(foreach dim,$(2),/mnt/cwd/$(1).K$(dim).fit.log)
endef
$(foreach sample,$(SAMPLES),$(eval $(call DIM_SUMMARY_RULE,$(sample),$(DIM_VALUES))))
$(foreach sample,$(NO_RIBO_MITO_SAMPLES),$(eval $(call DIM_SUMMARY_RULE,$(sample),$(NO_RIBO_MITO_DIM_VALUES))))

tls.dim_summary.pdf: tls1.dim_summary.tsv tls2.dim_summary.tsv \
		$(call DIM_FIT_GSEAS,tls1,$(DIM_VALUES)) $(call DIM_FIT_GSEAS,tls2,$(DIM_VALUES)) \
		$(SRC)/plot_combined_dim_summary.py $(SRC)/plot_dim_log_likelihood.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_combined_dim_summary.py \
		--sample TLS1 /mnt/cwd/tls1.dim_summary.tsv /mnt/cwd \
		--sample TLS2 /mnt/cwd/tls2.dim_summary.tsv /mnt/cwd \
		/mnt/cwd/$@

tls_no_ribo_mito.dim_summary.pdf: tls1_no_ribo_mito.dim_summary.tsv tls2_no_ribo_mito.dim_summary.tsv \
		$(call DIM_FIT_GSEAS,tls1_no_ribo_mito,$(NO_RIBO_MITO_DIM_VALUES)) $(call DIM_FIT_GSEAS,tls2_no_ribo_mito,$(NO_RIBO_MITO_DIM_VALUES)) \
		$(SRC)/plot_combined_dim_summary.py $(SRC)/plot_dim_log_likelihood.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_combined_dim_summary.py \
		--sample TLS1 /mnt/cwd/tls1_no_ribo_mito.dim_summary.tsv /mnt/cwd \
		--sample TLS2 /mnt/cwd/tls2_no_ribo_mito.dim_summary.tsv /mnt/cwd \
		/mnt/cwd/$@

clean-gexlineage:
	rm -f *.K*.fit.*
