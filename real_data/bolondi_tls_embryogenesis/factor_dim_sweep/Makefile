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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/factor_dim_sweep
PREPROCESSED_DIR := ../preprocessed_data
FILTERED_DIR := ../gex_filter

SAMPLES := tls1 tls2
DIM_VALUES ?= 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 30 40 50
DIM_FIT_SUMMARIES = $(foreach dim,$(DIM_VALUES),$(1).K$(dim).fit.summary.tsv)
DIM_FIT_GSEAS = $(foreach dim,$(DIM_VALUES),$(1).K$(dim).fit.L.gsea.tsv)
DIM_FIT_TIMES = $(foreach dim,$(DIM_VALUES),$(1).K$(dim).fit.time)

.PHONY: all clean-gexlineage
all: $(SAMPLES:%=%.dim_summary.pdf)

tls1.K%.fit.summary.tsv: $(FILTERED_DIR)/tls1.tree.nex $(FILTERED_DIR)/tls1.moran.filtered.tsv
	/usr/bin/time -o tls1.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls1.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls1.K$*.fit.term

tls2.K%.fit.summary.tsv: $(FILTERED_DIR)/tls2.tree.nex $(FILTERED_DIR)/tls2.moran.filtered.tsv
	/usr/bin/time -o tls2.K$*.fit.time $(GEXFACTOR) --trees $< --expr $(word 2,$^) --outprefix tls2.K$*.fit --dim $* --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > tls2.K$*.fit.term

tls1.K%.fit.L.tsv tls1.K%.fit.log tls1.K%.fit.time: tls1.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }
tls2.K%.fit.L.tsv tls2.K%.fit.log tls2.K%.fit.time: tls2.K%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

%.fit.L.gsea.tsv: %.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< mouse /mnt/cwd/$*.fit.L.gsea > $*.fit.gsea.term

tls1.dim_summary.pdf: $(call DIM_FIT_SUMMARIES,tls1) $(call DIM_FIT_GSEAS,tls1) $(call DIM_FIT_TIMES,tls1) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(DIM_VALUES),tls1.K$(dim).fit.log); do test -e $$log || { echo "ERROR: expected fit log missing: $$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_dim_log_likelihood.py /mnt/cwd/$@ $(foreach dim,$(DIM_VALUES),/mnt/cwd/tls1.K$(dim).fit.log)

tls2.dim_summary.pdf: $(call DIM_FIT_SUMMARIES,tls2) $(call DIM_FIT_GSEAS,tls2) $(call DIM_FIT_TIMES,tls2) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(foreach dim,$(DIM_VALUES),tls2.K$(dim).fit.log); do test -e $$log || { echo "ERROR: expected fit log missing: $$log" >&2; exit 1; }; done
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_dim_log_likelihood.py /mnt/cwd/$@ $(foreach dim,$(DIM_VALUES),/mnt/cwd/tls2.K$(dim).fit.log)

clean-gexlineage:
	rm -f *.K*.fit.*
