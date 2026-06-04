export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXLINEAGE := $(GEX_LINEAGE_DIR)/bin/gexLineage
BOLONDI_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis
PYTHON ?= /home/staklins/miniforge3/envs/pyplotting/bin/python

CLONES ?= tls1 tls2

# Ability to downsample the input data for testing. Leave these empty to use the full data.
DOWNSAMPLE_CELLS ?= 100
DOWNSAMPLE_GENES ?= 1000
DOWNSAMPLE_ENABLED := $(strip $(DOWNSAMPLE_CELLS)$(DOWNSAMPLE_GENES))
DOWNSAMPLE_CELLS_LABEL := $(if $(strip $(DOWNSAMPLE_CELLS)),$(DOWNSAMPLE_CELLS),all)
DOWNSAMPLE_GENES_LABEL := $(if $(strip $(DOWNSAMPLE_GENES)),$(DOWNSAMPLE_GENES),all)

# These lists define the sweep. Keep them short at first; the total number of
# runs is the product of clones and all value-list lengths.
L_L1_STRENGTHS ?= 0 5.0
L_CORRELATION_STRENGTHS ?= 0 5.0
F_CORRELATION_STRENGTHS ?= 0 5.0
F_ORTHOGONALITY_STRENGTHS ?= 0 5.0
FINAL_ABSORBING_VALUES ?= 0 1
FINAL_ABSORBING_L2_STRENGTHS ?= 1.0
VARIMAX_VALUES ?= 0 1

TREE_TOTAL_TIME ?= 1.0
DIM ?= 5
MAX_ITER ?= 10000
NTHREADS ?= 1
SEED ?= 1

BASE_GEXLINEAGE_OPTS := \
	--tree-total-time $(TREE_TOTAL_TIME) \
	--filter-test moran \
	--dim $(DIM) \
	--max-iter $(MAX_ITER) \
	--nthreads $(NTHREADS) \
	--verbose-log \
	--seed $(SEED)

define GEX_INPUT
$(if $(DOWNSAMPLE_ENABLED),runs/gex_inputs/$(1).cells$(DOWNSAMPLE_CELLS_LABEL).genes$(DOWNSAMPLE_GENES_LABEL).tsv,$(BOLONDI_DIR)/$(1)_gex.tsv)
endef

define CONDITION_NAME
Ll$(1)_Lc$(2)_Fc$(3)_Fo$(4)_Fa$(5)_Fl2$(6)_V$(7)
endef

define CONDITION_OPTS
--L-l1-strength $(1) \
$(if $(filter-out 0 0.0 0.00 0.000,$(2)),--L-correlation-strength $(2),) \
$(if $(filter-out 0 0.0 0.00 0.000,$(3)),--F-correlation-strength $(3),) \
$(if $(filter-out 0 0.0 0.00 0.000,$(4)),--F-orthogonality-strength $(4),) \
$(if $(filter 1 true yes,$(5)),--final-absorbing-factor --final-absorbing-L-l2-strength $(6),) \
$(if $(filter 1 true yes,$(7)),--varimax,)
endef

CONDITIONS := $(foreach ll1,$(L_L1_STRENGTHS),$(foreach lc,$(L_CORRELATION_STRENGTHS),$(foreach fc,$(F_CORRELATION_STRENGTHS),$(foreach fo,$(F_ORTHOGONALITY_STRENGTHS),$(foreach fa,$(FINAL_ABSORBING_VALUES),$(foreach fl2,$(FINAL_ABSORBING_L2_STRENGTHS),$(foreach vx,$(VARIMAX_VALUES),$(call CONDITION_NAME,$(ll1),$(lc),$(fc),$(fo),$(fa),$(fl2),$(vx)))))))))

RUN_PREFIXES := $(foreach clone,$(CLONES),$(foreach cond,$(CONDITIONS),runs/$(clone)/$(cond)))
RUN_SUMMARIES := $(addsuffix .summary.tsv,$(RUN_PREFIXES))
RUN_METRICS := $(addsuffix .overlap.tsv,$(RUN_PREFIXES))
PLOT_PDF := sweep_summary_plots.pdf

all: metrics
run: $(RUN_SUMMARIES)
metrics: overlap_summary.tsv
plot: $(PLOT_PDF)

print-config:
	@echo "CLONES=$(CLONES)"
	@echo "CONDITIONS=$(CONDITIONS)"
	@echo "L_L1_STRENGTHS=$(L_L1_STRENGTHS)"
	@echo "L_CORRELATION_STRENGTHS=$(L_CORRELATION_STRENGTHS)"
	@echo "F_CORRELATION_STRENGTHS=$(F_CORRELATION_STRENGTHS)"
	@echo "F_ORTHOGONALITY_STRENGTHS=$(F_ORTHOGONALITY_STRENGTHS)"
	@echo "FINAL_ABSORBING_VALUES=$(FINAL_ABSORBING_VALUES)"
	@echo "FINAL_ABSORBING_L2_STRENGTHS=$(FINAL_ABSORBING_L2_STRENGTHS)"
	@echo "VARIMAX_VALUES=$(VARIMAX_VALUES)"
	@echo "DOWNSAMPLE_CELLS=$(DOWNSAMPLE_CELLS)"
	@echo "DOWNSAMPLE_GENES=$(DOWNSAMPLE_GENES)"
	@echo "GEXLINEAGE=$(GEXLINEAGE)"
	@echo "BOLONDI_DIR=$(BOLONDI_DIR)"
	@echo "BASE_GEXLINEAGE_OPTS=$(BASE_GEXLINEAGE_OPTS)"

runs/gex_inputs/%.cells$(DOWNSAMPLE_CELLS_LABEL).genes$(DOWNSAMPLE_GENES_LABEL).tsv: $(BOLONDI_DIR)/%_gex.tsv
	mkdir -p runs/gex_inputs
	awk -F'\t' -v OFS='\t' -v max_cells="$(DOWNSAMPLE_CELLS)" -v max_genes="$(DOWNSAMPLE_GENES)" '\
		BEGIN { \
			if (max_cells != "" && max_cells !~ /^[0-9]+$$/) { \
				printf "DOWNSAMPLE_CELLS must be a non-negative integer, got %s\n", max_cells > "/dev/stderr"; \
				exit 1; \
			} \
			if (max_genes != "" && max_genes !~ /^[0-9]+$$/) { \
				printf "DOWNSAMPLE_GENES must be a non-negative integer, got %s\n", max_genes > "/dev/stderr"; \
				exit 1; \
			} \
		} \
		NR == 1 { \
			n_fields = NF; \
			if (max_genes != "" && max_genes + 1 < n_fields) n_fields = max_genes + 1; \
		} \
		NR > 1 && max_cells != "" && NR > max_cells + 1 { exit } \
		{ \
			for (i = 1; i <= n_fields; i++) { \
				printf "%s%s", (i == 1 ? "" : OFS), $$i; \
			} \
			printf "\n"; \
		}' $< > $@

define RUN_TEMPLATE
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).term \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).time \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).summary.tsv \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).F.tsv \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).L.tsv \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).X.tsv \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).log \
runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).L.pearson_correlation.tsv: \
$(BOLONDI_DIR)/$(1).vine.nex $(call GEX_INPUT,$(1))
	mkdir -p runs/$(1)
	/usr/bin/time -o runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).time $(GEXLINEAGE) \
		--trees $(BOLONDI_DIR)/$(1).vine.nex \
		--expr $(call GEX_INPUT,$(1)) \
		--outprefix runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)) \
		$(BASE_GEXLINEAGE_OPTS) \
		$(call CONDITION_OPTS,$(2),$(3),$(4),$(5),$(6),$(7),$(8)) > runs/$(1)/$(call CONDITION_NAME,$(2),$(3),$(4),$(5),$(6),$(7),$(8)).term
endef

$(foreach clone,$(CLONES),$(foreach ll1,$(L_L1_STRENGTHS),$(foreach lc,$(L_CORRELATION_STRENGTHS),$(foreach fc,$(F_CORRELATION_STRENGTHS),$(foreach fo,$(F_ORTHOGONALITY_STRENGTHS),$(foreach fa,$(FINAL_ABSORBING_VALUES),$(foreach fl2,$(FINAL_ABSORBING_L2_STRENGTHS),$(foreach vx,$(VARIMAX_VALUES),$(eval $(call RUN_TEMPLATE,$(clone),$(ll1),$(lc),$(fc),$(fo),$(fa),$(fl2),$(vx)))))))))))

%.overlap.tsv: %.L.pearson_correlation.tsv %.summary.tsv %.log
	prefix="$*"; \
	clone=$$(basename $$(dirname "$$prefix")); \
	condition=$$(basename "$$prefix"); \
	corr_metrics=$$(awk -F'\t' '\
		NR == 1 { next } \
		{ \
			row = NR - 1; \
			for (i = 2; i <= NF; i++) { \
				col = i - 1; \
				if (row != col) { \
					v = $$i + 0; \
					abs_v = v; \
					if (abs_v < 0) abs_v = -abs_v; \
					sum_abs += abs_v; \
					sum += v; \
					n++; \
				} \
			} \
		} \
		END { \
			if (n > 0) printf "%.17g\t%.17g\t%d", sum_abs / n, sum / n, n; \
			else printf "NA\tNA\t0"; \
		}' $<); \
	summary_metrics=$$(awk -F'\t' '\
		$$1 == "brownian_negll" { brownian = $$2 } \
		$$1 == "observation_negll" { observation = $$2 } \
		END { \
			if (brownian == "") brownian = "NA"; \
			if (observation == "") observation = "NA"; \
			printf "%s\t%s", brownian, observation; \
		}' "$$prefix.summary.tsv"); \
	final_objective=$$(awk -F'\t' '$$1 ~ /^[0-9]+$$/ { objective = $$2 } END { if (objective == "") objective = "NA"; print objective }' "$$prefix.log"); \
	printf "%s\t%s\t%s\t%s\t%s\n" "$$clone" "$$condition" "$$corr_metrics" "$$final_objective" "$$summary_metrics" > $@

overlap_summary.tsv: $(RUN_METRICS)
	echo -e "clone\tcondition\tmean_offdiag_abs_L_correlation\tmean_offdiag_L_correlation\tn_offdiag_entries\tfinal_objective\tbrownian_negll\tobservation_negll" > $@
	cat $(RUN_METRICS) >> $@

$(PLOT_PDF): overlap_summary.tsv plot_sweep_summary.py
	$(PYTHON) plot_sweep_summary.py overlap_summary.tsv $@

clean-metrics:
	rm -f overlap_summary.tsv $(RUN_METRICS)

clean-plots:
	rm -f $(PLOT_PDF)

clean:
	rm -rf runs overlap_summary.tsv $(PLOT_PDF)
