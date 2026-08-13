export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:
.SECONDEXPANSION:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
GEXFACTOR := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexFactor
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_dim_sweep
PREPROCESSED_DIR := ../preprocessed_data
FILTERED_DIR := ../gex_filter

CLONES := 1 2 3 4
DIM_VALUES ?= 1 2 3 4 5 10 15
NO_RIBO_MITO_DIM_VALUES ?= $(DIM_VALUES) 20 25 30
FIT_CP = $(firstword $(subst .K, ,$(1)))
FIT_DIM = $(word 2,$(subst .K, ,$(1)))
FIT_TREE = $(PREPROCESSED_DIR)/quinn_$(call FIT_CP,$(1)).vine.mean.nex
FIT_GEX = $(FILTERED_DIR)/quinn_$(call FIT_CP,$(1)).moran.filtered.tsv
DIM_FIT_SUMMARIES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.summary.tsv)
DIM_FIT_GSEAS = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.L.gsea.tsv)
DIM_FIT_TIMES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.time)
NO_RIBO_MITO_DIM_FIT_SUMMARIES = $(foreach dim,$(NO_RIBO_MITO_DIM_VALUES),quinn_$(1)_no_ribo_mito.K$(dim).fit.summary.tsv)
NO_RIBO_MITO_DIM_FIT_GSEAS = $(foreach dim,$(NO_RIBO_MITO_DIM_VALUES),quinn_$(1)_no_ribo_mito.K$(dim).fit.L.gsea.tsv)
NO_RIBO_MITO_DIM_FIT_TIMES = $(foreach dim,$(NO_RIBO_MITO_DIM_VALUES),quinn_$(1)_no_ribo_mito.K$(dim).fit.time)

.PHONY: all archive clean-gexlineage
all: quinn.dim_summary.pdf quinn_no_ribo_mito.dim_summary.pdf

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

quinn_%.fit.summary.tsv: $$(call FIT_TREE,$$*) $$(call FIT_GEX,$$*)
	/usr/bin/time -o quinn_$*.fit.time $(GEXFACTOR) --trees $(call FIT_TREE,$*) --expr $(call FIT_GEX,$*) --outprefix quinn_$*.fit --dim $(call FIT_DIM,$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_$*.fit.term

quinn_1_no_ribo_mito.K%.fit.summary.tsv: $(PREPROCESSED_DIR)/quinn_1.vine.mean.nex $(FILTERED_DIR)/quinn_1_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o quinn_1_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix quinn_1_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_1_no_ribo_mito.K$*.fit.term

quinn_2_no_ribo_mito.K%.fit.summary.tsv: $(PREPROCESSED_DIR)/quinn_2.vine.mean.nex $(FILTERED_DIR)/quinn_2_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o quinn_2_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix quinn_2_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_2_no_ribo_mito.K$*.fit.term

quinn_3_no_ribo_mito.K%.fit.summary.tsv: $(PREPROCESSED_DIR)/quinn_3.vine.mean.nex $(FILTERED_DIR)/quinn_3_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o quinn_3_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix quinn_3_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_3_no_ribo_mito.K$*.fit.term

quinn_4_no_ribo_mito.K%.fit.summary.tsv: $(PREPROCESSED_DIR)/quinn_4.vine.mean.nex $(FILTERED_DIR)/quinn_4_no_ribo_mito.moran.filtered.tsv
	/usr/bin/time -o quinn_4_no_ribo_mito.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix quinn_4_no_ribo_mito.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_4_no_ribo_mito.K$*.fit.term

quinn_%.fit.L.tsv quinn_%.fit.log quinn_%.fit.time: quinn_%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

quinn_%.fit.L.gsea.tsv: quinn_%.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< human /mnt/cwd/quinn_$*.fit.L.gsea > quinn_$*.fit.gsea.term

define DIM_SUMMARY_RULE
quinn_$(1).dim_summary.tsv: $(call DIM_FIT_SUMMARIES,$(1)) $(call DIM_FIT_GSEAS,$(1)) \
		$(call DIM_FIT_TIMES,$(1)) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.log); do test -e $$$$log || { echo "ERROR: expected fit log missing: $$$$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_dim_log_likelihood.py --summary-only \
		/mnt/cwd/quinn_$(1).dim_summary.pdf \
		$(foreach dim,$(DIM_VALUES),/mnt/cwd/quinn_$(1).K$(dim).fit.log)
endef
$(foreach cp,$(CLONES),$(eval $(call DIM_SUMMARY_RULE,$(cp))))

define NO_RIBO_MITO_DIM_SUMMARY_RULE
quinn_$(1)_no_ribo_mito.dim_summary.tsv: $(call NO_RIBO_MITO_DIM_FIT_SUMMARIES,$(1)) $(call NO_RIBO_MITO_DIM_FIT_GSEAS,$(1)) \
		$(call NO_RIBO_MITO_DIM_FIT_TIMES,$(1)) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(NO_RIBO_MITO_DIM_VALUES),quinn_$(1)_no_ribo_mito.K$(dim).fit.log); do test -e $$$$log || { echo "ERROR: expected fit log missing: $$$$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_dim_log_likelihood.py --summary-only \
		/mnt/cwd/quinn_$(1)_no_ribo_mito.dim_summary.pdf \
		$(foreach dim,$(NO_RIBO_MITO_DIM_VALUES),/mnt/cwd/quinn_$(1)_no_ribo_mito.K$(dim).fit.log)
endef
$(foreach cp,$(CLONES),$(eval $(call NO_RIBO_MITO_DIM_SUMMARY_RULE,$(cp))))

quinn.dim_summary.pdf: $(CLONES:%=quinn_%.dim_summary.tsv) \
		$(SRC)/plot_combined_dim_summary.py $(SRC)/plot_dim_log_likelihood.py
	@for gsea in $(foreach cp,$(CLONES),$(call DIM_FIT_GSEAS,$(cp))); do test -e $$gsea || { echo "ERROR: existing GSEA result missing: $$gsea" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_combined_dim_summary.py \
		$(foreach cp,$(CLONES),--sample CP$(cp) /mnt/cwd/quinn_$(cp).dim_summary.tsv /mnt/cwd) \
		/mnt/cwd/$@

quinn_no_ribo_mito.dim_summary.pdf: $(CLONES:%=quinn_%_no_ribo_mito.dim_summary.tsv) \
		$(SRC)/plot_combined_dim_summary.py $(SRC)/plot_dim_log_likelihood.py
	@for gsea in $(foreach cp,$(CLONES),$(call NO_RIBO_MITO_DIM_FIT_GSEAS,$(cp))); do test -e $$gsea || { echo "ERROR: existing GSEA result missing: $$gsea" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_combined_dim_summary.py \
		$(foreach cp,$(CLONES),--sample CP$(cp) /mnt/cwd/quinn_$(cp)_no_ribo_mito.dim_summary.tsv /mnt/cwd) \
		/mnt/cwd/$@

clean-gexlineage:
	rm -f *.K*.fit.*
