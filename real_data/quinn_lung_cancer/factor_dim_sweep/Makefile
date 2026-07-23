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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/factor_dim_sweep
PREPROCESSED_DIR := ../preprocessed_data
FILTERED_DIR := ../gex_filter

CLONES := 1 2 3 4
DIM_VALUES ?= 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 30 40 50
FIT_CP = $(firstword $(subst .K, ,$(1)))
FIT_DIM = $(word 2,$(subst .K, ,$(1)))
FIT_TREE = $(FILTERED_DIR)/quinn_$(call FIT_CP,$(1)).tree.nex
FIT_GEX = $(FILTERED_DIR)/quinn_$(call FIT_CP,$(1)).moran.filtered.tsv
DIM_FIT_SUMMARIES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.summary.tsv)
DIM_FIT_GSEAS = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.L.gsea.tsv)
DIM_FIT_TIMES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.time)

.PHONY: all clean-gexlineage
all: $(CLONES:%=quinn_%.dim_summary.pdf)

quinn_%.fit.summary.tsv: $$(call FIT_TREE,$$*) $$(call FIT_GEX,$$*)
	/usr/bin/time -o quinn_$*.fit.time $(GEXFACTOR) --trees $(call FIT_TREE,$*) --expr $(call FIT_GEX,$*) --outprefix quinn_$*.fit --dim $(call FIT_DIM,$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_$*.fit.term

quinn_%.fit.L.tsv quinn_%.fit.log quinn_%.fit.time: quinn_%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

quinn_%.fit.L.gsea.tsv: quinn_%.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< human /mnt/cwd/quinn_$*.fit.L.gsea > quinn_$*.fit.gsea.term

quinn_%.dim_summary.pdf: $$(call DIM_FIT_SUMMARIES,$$*) $$(call DIM_FIT_GSEAS,$$*) $$(call DIM_FIT_TIMES,$$*) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(DIM_VALUES),quinn_$*.K$(dim).fit.log); do test -e $$log || { echo "ERROR: expected fit log missing: $$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_dim_log_likelihood.py /mnt/cwd/$@ $(foreach dim,$(DIM_VALUES),/mnt/cwd/quinn_$*.K$(dim).fit.log)

clean-gexlineage:
	rm -f *.K*.fit.*
