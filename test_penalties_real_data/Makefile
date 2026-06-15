export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXLINEAGE := $(GEX_LINEAGE_DIR)/bin/gexLineage
BOLONDI_DIR := $(BENCHMARKS_DIR)/bolondi_tls_embryogenesis
QUINN_DIR := $(BENCHMARKS_DIR)/quinn_lung_cancer
PYTHON ?= /home/staklins/miniforge3/envs/pyplotting/bin/python

# CLONES can include TLS clones, such as tls1/tls2, and Quinn shorthand clones,
# such as quinn30 for Quinn CP 30.
# CLONES ?= tls1 tls2
CLONES ?= quinn10

# Ability to downsample the input data for testing. Leave these empty to use the full data.
DOWNSAMPLE_CELLS ?= 
DOWNSAMPLE_GENES ?= 
DOWNSAMPLE_ENABLED := $(strip $(DOWNSAMPLE_CELLS)$(DOWNSAMPLE_GENES))
DOWNSAMPLE_CELLS_LABEL := $(if $(strip $(DOWNSAMPLE_CELLS)),$(DOWNSAMPLE_CELLS),all)
DOWNSAMPLE_GENES_LABEL := $(if $(strip $(DOWNSAMPLE_GENES)),$(DOWNSAMPLE_GENES),all)

# These lists define the sweep. Keep them short at first; the total number of
# runs is the product of clones and all value-list lengths.
L_L1_STRENGTHS ?= 5.0
L_CORRELATION_STRENGTHS ?= 0
L_LOADING_OVERLAP_STRENGTHS ?= 100
F_CORRELATION_STRENGTHS ?= 0
F_ORTHOGONALITY_STRENGTHS ?= 0
FINAL_ABSORBING_VALUES ?= none
FINAL_ABSORBING_L2_STRENGTHS ?= 0
VARIMAX_VALUES ?= 0

TREE_MODES ?= one
TREE_INDEX ?= 1
RIBO_MITO_VALUES ?= 1
TREE_TOTAL_TIME ?= 1.0
DIM ?= 5
DIM_VALUES ?= $(DIM)
MAX_ITER ?= 20000
NTHREADS ?= 1
SEED ?= 1

BASE_GEXLINEAGE_OPTS := \
	--tree-total-time $(TREE_TOTAL_TIME) \
	--filter-test moran \
	--max-iter $(MAX_ITER) \
	--nthreads $(NTHREADS) \
	--verbose-log \
	--seed $(SEED)

define IS_QUINN_CLONE
$(filter quinn%,$(1))
endef

define QUINN_CP
$(patsubst quinn%,%,$(1))
endef

define RAW_GEX_INPUT
$(if $(call IS_QUINN_CLONE,$(1)),$(QUINN_DIR)/quinn_$(call QUINN_CP,$(1))_gex.tsv,$(BOLONDI_DIR)/$(1)_gex.tsv)
endef

define TREE_INPUT
$(if $(call IS_QUINN_CLONE,$(1)),$(QUINN_DIR)/quinn_$(call QUINN_CP,$(1)).vine.nex,$(BOLONDI_DIR)/$(1).vine.nex)
endef

define TREE_MODE_NAME
$(if $(filter one,$(1)),one$(TREE_INDEX),$(1))
endef

define RUN_TREE_INPUT
$(call TREE_INPUT,$(1))
endef

define TREE_MODE_OPTS
$(if $(filter one,$(1)),--tree-index $(TREE_INDEX),)
endef

define RIBO_MITO_OPTS
$(if $(filter 1 true yes,$(1)),--remove-ribo-mito-genes,)
endef

define GEX_INPUT
$(if $(DOWNSAMPLE_ENABLED),runs/gex_inputs/$(1).cells$(DOWNSAMPLE_CELLS_LABEL).genes$(DOWNSAMPLE_GENES_LABEL).tsv,$(call RAW_GEX_INPUT,$(1)))
endef

define CONDITION_NAME
T$(1)_Rm$(2)_K$(3)_Ll$(4)_Lc$(5)_Lo$(6)_Fc$(7)_Fo$(8)_Fa$(9)_Fl2$(10)_V$(11)
endef

define CONDITION_OPTS
--L-l1-strength $(1) \
$(if $(filter-out 0 0.0 0.00 0.000,$(2)),--L-correlation-strength $(2),) \
$(if $(filter-out 0 0.0 0.00 0.000,$(3)),--L-loading-overlap-strength $(3),) \
$(if $(filter-out 0 0.0 0.00 0.000,$(4)),--F-correlation-strength $(4),) \
$(if $(filter-out 0 0.0 0.00 0.000,$(5)),--F-orthogonality-strength $(5),) \
$(if $(filter-out none 0 0.0 false no,$(6)),--final-absorbing-factor $(6) --final-absorbing-L-l2-strength $(7),) \
$(if $(filter 1 true yes,$(8)),--varimax,)
endef

define ABSORBING_L2_VALUES
$(if $(filter none 0 0.0 false no,$(1)),0,$(FINAL_ABSORBING_L2_STRENGTHS))
endef

CONDITIONS := $(foreach tm,$(TREE_MODES),$(foreach rm,$(RIBO_MITO_VALUES),$(foreach k,$(DIM_VALUES),$(foreach ll1,$(L_L1_STRENGTHS),$(foreach lc,$(L_CORRELATION_STRENGTHS),$(foreach lo,$(L_LOADING_OVERLAP_STRENGTHS),$(foreach fc,$(F_CORRELATION_STRENGTHS),$(foreach fo,$(F_ORTHOGONALITY_STRENGTHS),$(foreach fa,$(FINAL_ABSORBING_VALUES),$(foreach fl2,$(call ABSORBING_L2_VALUES,$(fa)),$(foreach vx,$(VARIMAX_VALUES),$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(tm)),$(rm),$(k),$(ll1),$(lc),$(lo),$(fc),$(fo),$(fa),$(fl2),$(vx)))))))))))))

define CLONE_RUN_PREFIXES
$(foreach cond,$(CONDITIONS),runs/$(1)/$(cond))
endef

RUN_PREFIXES := $(foreach clone,$(CLONES),$(call CLONE_RUN_PREFIXES,$(clone)))
RUN_SUMMARIES := $(addsuffix .summary.tsv,$(RUN_PREFIXES))
RUN_METRICS := $(addsuffix .overlap.tsv,$(RUN_PREFIXES))
RUN_TIMES := $(addsuffix .time,$(RUN_PREFIXES))
RUN_LOG_PDFS := $(addsuffix .log.pdf,$(RUN_PREFIXES))
RUN_FACTOR_HEATMAPS := $(addsuffix .L.factor_heatmap.pdf,$(RUN_PREFIXES))
RUN_FLOWTREES := $(addsuffix .latent_flow.factor_trees.pdf,$(RUN_PREFIXES))
RUN_FLOW_DISTANCES := $(addsuffix .latent_flow.total_euclidean_distance.tsv,$(RUN_PREFIXES))
RUN_GSEAS := $(addsuffix .L.gsea.tsv,$(RUN_PREFIXES))
RUN_GSEA_TERMS := $(addsuffix .gsea.term,$(RUN_PREFIXES))
RUN_TOP_GENES := $(addsuffix .L.top_genes.tsv,$(RUN_PREFIXES))
SUMMARY_DIRS := $(foreach clone,$(CLONES),runs/$(clone)_summary)
L_CORR_SUMMARY_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/summary.l_correlation.$(clone).tsv)
GSEA_SUMMARY_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/summary.gsea_jaccard.$(clone).tsv)
TOP_GENE_SUMMARY_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/summary.top_genes_jaccard.$(clone).tsv)
RUNTIME_SUMMARY_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/summary.runtime.$(clone).tsv)
L_CORR_PLOT_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/sweep_summary.l_correlation.$(clone).pdf)
GSEA_PLOT_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/sweep_summary.gsea_jaccard.$(clone).pdf)
TOP_GENE_PLOT_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/sweep_summary.top_genes_jaccard.$(clone).pdf)
RUNTIME_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/runtime_marginals.$(clone).pdf)
PENALTY_SELECTION_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/selected_penalty_runs.l_correlation.$(clone).tsv)
PENALTY_WORST_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/worst_penalty_runs.l_correlation.$(clone).tsv)
PENALTY_SELECTION_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/penalty_selection.l_correlation.$(clone).pdf)
PENALTY_GSEA_SELECTION_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/selected_penalty_runs.gsea_jaccard.$(clone).tsv)
PENALTY_GSEA_WORST_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/worst_penalty_runs.gsea_jaccard.$(clone).tsv)
PENALTY_GSEA_SELECTION_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/penalty_selection.gsea_jaccard.$(clone).pdf)
PENALTY_TOP_GENE_SELECTION_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/selected_penalty_runs.top_genes_jaccard.$(clone).tsv)
PENALTY_TOP_GENE_WORST_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/worst_penalty_runs.top_genes_jaccard.$(clone).tsv)
PENALTY_TOP_GENE_SELECTION_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/penalty_selection.top_genes_jaccard.$(clone).pdf)
OLD_SUMMARY_OUTPUTS := overlap_summary.tsv overlap_summary_gsea.tsv overlap_summary_top_genes.tsv penalty_runtime_summary.tsv
OLD_PLOT_OUTPUTS := sweep_summary_plots.pdf sweep_summary_gsea_plots.pdf sweep_summary_top_gene_plots.pdf penalty_runtime_marginals.pdf penalty_selection_plots.pdf penalty_selection_gsea_plots.pdf penalty_selection_top_gene_plots.pdf
OLD_CLONE_OUTPUTS := $(foreach clone,$(CLONES),overlap_summary.$(clone).tsv overlap_summary_gsea.$(clone).tsv overlap_summary_top_genes.$(clone).tsv summary.l_correlation.$(clone).tsv summary.gsea.$(clone).tsv summary.gsea_jaccard.$(clone).tsv summary.top_genes.$(clone).tsv summary.top_genes_jaccard.$(clone).tsv summary.runtime.$(clone).tsv sweep_summary_plots.$(clone).pdf sweep_summary_gsea_plots.$(clone).pdf sweep_summary_top_gene_plots.$(clone).pdf sweep_summary.l_correlation.$(clone).pdf sweep_summary.gsea.$(clone).pdf sweep_summary.gsea_jaccard.$(clone).pdf sweep_summary.top_genes.$(clone).pdf sweep_summary.top_genes_jaccard.$(clone).pdf runtime_marginals.$(clone).pdf selected_penalty_runs.$(clone).tsv worst_penalty_runs.$(clone).tsv penalty_selection_plots.$(clone).pdf selected_penalty_gsea_runs.$(clone).tsv worst_penalty_gsea_runs.$(clone).tsv penalty_selection_gsea_plots.$(clone).pdf selected_penalty_top_gene_runs.$(clone).tsv worst_penalty_top_gene_runs.$(clone).tsv penalty_selection_top_gene_plots.$(clone).pdf selected_penalty_runs.l_correlation.$(clone).tsv worst_penalty_runs.l_correlation.$(clone).tsv penalty_selection.l_correlation.$(clone).pdf selected_penalty_runs.gsea.$(clone).tsv worst_penalty_runs.gsea.$(clone).tsv penalty_selection.gsea.$(clone).pdf selected_penalty_runs.gsea_jaccard.$(clone).tsv worst_penalty_runs.gsea_jaccard.$(clone).tsv penalty_selection.gsea_jaccard.$(clone).pdf selected_penalty_runs.top_genes.$(clone).tsv worst_penalty_runs.top_genes.$(clone).tsv penalty_selection.top_genes.$(clone).pdf selected_penalty_runs.top_genes_jaccard.$(clone).tsv worst_penalty_runs.top_genes_jaccard.$(clone).tsv penalty_selection.top_genes_jaccard.$(clone).pdf)
GSEA_PADJ_CUTOFF ?= 1.0
GSEA_TOP_N ?= 20

all: run metrics plots
run: $(RUN_SUMMARIES)
metrics: $(L_CORR_SUMMARY_TSVS) $(GSEA_SUMMARY_TSVS) $(TOP_GENE_SUMMARY_TSVS) $(RUNTIME_SUMMARY_TSVS) $(PENALTY_SELECTION_TSVS) $(PENALTY_WORST_TSVS) $(RUN_GSEAS)
plots: $(RUN_LOG_PDFS) $(RUN_FACTOR_HEATMAPS) $(RUN_FLOWTREES) $(L_CORR_PLOT_PDFS) $(GSEA_PLOT_PDFS) $(TOP_GENE_PLOT_PDFS) $(RUNTIME_PDFS) $(PENALTY_SELECTION_PDFS) $(PENALTY_GSEA_SELECTION_PDFS) $(PENALTY_TOP_GENE_SELECTION_PDFS)

print-config:
	@echo "CLONES=$(CLONES)"
	@echo "CONDITIONS=$(CONDITIONS)"
	@echo "L_L1_STRENGTHS=$(L_L1_STRENGTHS)"
	@echo "L_CORRELATION_STRENGTHS=$(L_CORRELATION_STRENGTHS)"
	@echo "L_LOADING_OVERLAP_STRENGTHS=$(L_LOADING_OVERLAP_STRENGTHS)"
	@echo "F_CORRELATION_STRENGTHS=$(F_CORRELATION_STRENGTHS)"
	@echo "F_ORTHOGONALITY_STRENGTHS=$(F_ORTHOGONALITY_STRENGTHS)"
	@echo "FINAL_ABSORBING_VALUES=$(FINAL_ABSORBING_VALUES)"
	@echo "FINAL_ABSORBING_L2_STRENGTHS=$(FINAL_ABSORBING_L2_STRENGTHS)"
	@echo "VARIMAX_VALUES=$(VARIMAX_VALUES)"
	@echo "TREE_MODES=$(TREE_MODES)"
	@echo "TREE_INDEX=$(TREE_INDEX)"
	@echo "RIBO_MITO_VALUES=$(RIBO_MITO_VALUES)"
	@echo "DIM_VALUES=$(DIM_VALUES)"
	@echo "DOWNSAMPLE_CELLS=$(DOWNSAMPLE_CELLS)"
	@echo "DOWNSAMPLE_GENES=$(DOWNSAMPLE_GENES)"
	@echo "GSEA_PADJ_CUTOFF=$(GSEA_PADJ_CUTOFF)"
	@echo "GSEA_TOP_N=$(GSEA_TOP_N)"
	@echo "GEXLINEAGE=$(GEXLINEAGE)"
	@echo "BOLONDI_DIR=$(BOLONDI_DIR)"
	@echo "QUINN_DIR=$(QUINN_DIR)"
	@echo "PYPLOTTING_SIF=$(PYPLOTTING_SIF)"
	@echo "BASE_GEXLINEAGE_OPTS=$(BASE_GEXLINEAGE_OPTS)"

define GEX_INPUT_TEMPLATE
runs/gex_inputs/$(1).cells$(DOWNSAMPLE_CELLS_LABEL).genes$(DOWNSAMPLE_GENES_LABEL).tsv: $(call RAW_GEX_INPUT,$(1))
	mkdir -p runs/gex_inputs
	awk -F'\t' -v OFS='\t' -v max_cells="$(DOWNSAMPLE_CELLS)" -v max_genes="$(DOWNSAMPLE_GENES)" '\
		BEGIN { \
			if (max_cells != "" && max_cells !~ /^[0-9]+$$$$/) { \
				printf "DOWNSAMPLE_CELLS must be a non-negative integer, got %s\n", max_cells > "/dev/stderr"; \
				exit 1; \
			} \
			if (max_genes != "" && max_genes !~ /^[0-9]+$$$$/) { \
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
				printf "%s%s", (i == 1 ? "" : OFS), $$$$i; \
			} \
			printf "\n"; \
		}' $$< > $$@
endef

$(foreach clone,$(CLONES),$(eval $(call GEX_INPUT_TEMPLATE,$(clone))))

define RUN_TEMPLATE
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).summary.tsv: \
$(call RUN_TREE_INPUT,$(1),$(2)) $(call GEX_INPUT,$(1))
	mkdir -p runs/$(1)
	/usr/bin/time -o runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).time $(GEXLINEAGE) \
		--trees $(call RUN_TREE_INPUT,$(1),$(2)) \
		--expr $(call GEX_INPUT,$(1)) \
		--outprefix runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)) \
		$(BASE_GEXLINEAGE_OPTS) \
		--dim $(4) \
		$(call TREE_MODE_OPTS,$(2)) \
		$(call RIBO_MITO_OPTS,$(3)) \
		$(call CONDITION_OPTS,$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)) > runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).term

runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).term \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).time \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).F.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).L.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).L.top_genes.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).X.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).log \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).latent_flow.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).L.pearson_correlation.tsv: \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6),$(7),$(8),$(9),$(10),$(11),$(12)).summary.tsv
	@test -e $$@ || { echo "ERROR: expected output missing after gexLineage run: $$@" >&2; exit 1; }
endef

$(foreach clone,$(CLONES),$(foreach tm,$(TREE_MODES),$(foreach rm,$(RIBO_MITO_VALUES),$(foreach k,$(DIM_VALUES),$(foreach ll1,$(L_L1_STRENGTHS),$(foreach lc,$(L_CORRELATION_STRENGTHS),$(foreach lo,$(L_LOADING_OVERLAP_STRENGTHS),$(foreach fc,$(F_CORRELATION_STRENGTHS),$(foreach fo,$(F_ORTHOGONALITY_STRENGTHS),$(foreach fa,$(FINAL_ABSORBING_VALUES),$(foreach fl2,$(call ABSORBING_L2_VALUES,$(fa)),$(foreach vx,$(VARIMAX_VALUES),$(eval $(call RUN_TEMPLATE,$(clone),$(tm),$(rm),$(k),$(ll1),$(lc),$(lo),$(fc),$(fo),$(fa),$(fl2),$(vx)))))))))))))))

define RUNTIME_TEMPLATE
runs/$(1)_summary/summary.runtime.$(1).tsv: $(foreach cond,$(CONDITIONS),runs/$(1)/$(cond).time) plot_penalty_runtime_marginals.py
	mkdir -p runs/$(1)_summary
	$(PYTHON) plot_penalty_runtime_marginals.py runs/$(1)_summary/summary.runtime.$(1).tsv runs/$(1)_summary/runtime_marginals.$(1).pdf $(foreach cond,$(CONDITIONS),runs/$(1)/$(cond).time)

runs/$(1)_summary/runtime_marginals.$(1).pdf: runs/$(1)_summary/summary.runtime.$(1).tsv
	@test -e $$@ || { echo "ERROR: expected output missing: $$@" >&2; exit 1; }
endef

$(foreach clone,$(CLONES),$(eval $(call RUNTIME_TEMPLATE,$(clone))))

%.overlap.tsv: %.L.pearson_correlation.tsv %.log
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
	fit_metrics=$$(awk -F'\t' '\
		BEGIN { objective = "NA"; brownian = "NA"; observation = "NA"; best_objective = ""; source = "NA"; } \
		$$1 == "step" { \
			for (i = 1; i <= NF; i++) col[$$i] = i; \
			next; \
		} \
		$$1 == "# best_state" && ("objective" in col) { \
			objective = $$(col["objective"] + 1); \
			if ("brownian_neglprior" in col) brownian = $$(col["brownian_neglprior"] + 1); \
			if ("observation_negll" in col) observation = $$(col["observation_negll"] + 1); \
			source = "best_state"; \
			next; \
		} \
		$$1 ~ /^[0-9]+$$/ && ("objective" in col) { \
			row_objective = $$(col["objective"]) + 0; \
			if (best_objective == "" || row_objective < best_objective) { \
				best_objective = row_objective; \
				objective = $$(col["objective"]); \
				if ("brownian_neglprior" in col) brownian = $$(col["brownian_neglprior"]); \
				if ("observation_negll" in col) observation = $$(col["observation_negll"]); \
				source = "best_sample"; \
			} \
		} \
		END { printf "%s\t%s\t%s\t%s", objective, brownian, observation, source; }' "$$prefix.log"); \
	printf "%s\t%s\t%s\t%s\n" "$$clone" "$$condition" "$$corr_metrics" "$$fit_metrics" > $@

%.log.pdf: %.log
	$(PYTHON) $(SRC)/plot_optimization_log.py $< -1 0 $<

%.L.factor_heatmap.pdf: %.L.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(BENCHMARKS_DIR):/mnt/benchmarks $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_heatmap.py \
			/mnt/benchmarks/test_penalties_real_data/$*.L.tsv \
			/mnt/benchmarks/test_penalties_real_data/$*.L.factor_heatmap.pdf

%.latent_flow.factor_trees.pdf %.latent_flow.total_euclidean_distance.tsv: %.latent_flow.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(BENCHMARKS_DIR):/mnt/benchmarks $(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow.py \
			/mnt/benchmarks/test_penalties_real_data/$*.latent_flow.tsv \
			/mnt/benchmarks/test_penalties_real_data/$*.latent_flow

%.L.gsea.tsv: %.L.tsv
	clone=$$(basename $$(dirname "$*")); \
	organism=mouse; \
	if [[ "$$clone" == quinn* ]]; then organism=human; fi; \
	singularity run --bind $(BENCHMARKS_DIR):/mnt/benchmarks $(CLUSTERPROFILER_SIF) \
		/mnt/benchmarks/src/run_factor_gsea.R \
		/mnt/benchmarks/test_penalties_real_data/$*.L.tsv \
		$$organism \
		/mnt/benchmarks/test_penalties_real_data/$*.L.gsea > $*.gsea.term

define SUMMARY_TEMPLATE
runs/$(1)_summary/summary.l_correlation.$(1).tsv: $(addsuffix .overlap.tsv,$(call CLONE_RUN_PREFIXES,$(1))) summarize_overlap_metrics.py
	mkdir -p runs/$(1)_summary
	$(PYTHON) summarize_overlap_metrics.py --runs-dir runs --output $$@ --clone $(1)

runs/$(1)_summary/sweep_summary.l_correlation.$(1).pdf: runs/$(1)_summary/summary.l_correlation.$(1).tsv plot_sweep_summary.py
	$(PYTHON) plot_sweep_summary.py $$< $$@ --clone $(1)

runs/$(1)_summary/sweep_summary.gsea_jaccard.$(1).pdf: runs/$(1)_summary/summary.gsea_jaccard.$(1).tsv plot_sweep_summary.py
	$(PYTHON) plot_sweep_summary.py $$< $$@ \
		--clone $(1) \
		--metric-col gsea_uniqueness \
		--metric-label "GSEA uniqueness (1 - mean Jaccard)"

runs/$(1)_summary/sweep_summary.top_genes_jaccard.$(1).pdf: runs/$(1)_summary/summary.top_genes_jaccard.$(1).tsv plot_sweep_summary.py
	$(PYTHON) plot_sweep_summary.py $$< $$@ \
		--clone $(1) \
		--metric-col top_gene_mean_pairwise_jaccard \
		--metric-label "Mean pairwise Jaccard of top up/down genes"

runs/$(1)_summary/selected_penalty_runs.l_correlation.$(1).tsv: runs/$(1)_summary/summary.l_correlation.$(1).tsv summarize_penalty_selection.py
	$(PYTHON) summarize_penalty_selection.py $$< runs/$(1)_summary/selected_penalty_runs.l_correlation.$(1).tsv runs/$(1)_summary/penalty_selection.l_correlation.$(1).pdf \
		--worst-tsv runs/$(1)_summary/worst_penalty_runs.l_correlation.$(1).tsv \
		--clone $(1)

runs/$(1)_summary/worst_penalty_runs.l_correlation.$(1).tsv runs/$(1)_summary/penalty_selection.l_correlation.$(1).pdf: runs/$(1)_summary/selected_penalty_runs.l_correlation.$(1).tsv
	@test -e $$@ || { echo "ERROR: expected output missing: $$@" >&2; exit 1; }

runs/$(1)_summary/summary.gsea_jaccard.$(1).tsv: runs/$(1)_summary/summary.l_correlation.$(1).tsv $(addsuffix .L.gsea.tsv,$(call CLONE_RUN_PREFIXES,$(1))) summarize_gsea_uniqueness.py
	$(PYTHON) summarize_gsea_uniqueness.py --runs-dir runs --base-summary $$< --output $$@ \
		--clone $(1) \
		--padj-cutoff $(GSEA_PADJ_CUTOFF) \
		--top-n $(GSEA_TOP_N)

runs/$(1)_summary/summary.top_genes_jaccard.$(1).tsv: runs/$(1)_summary/summary.l_correlation.$(1).tsv $(addsuffix .L.top_genes.tsv,$(call CLONE_RUN_PREFIXES,$(1))) summarize_top_gene_jaccard.py
	$(PYTHON) summarize_top_gene_jaccard.py --runs-dir runs --base-summary $$< --output $$@ --clone $(1)

runs/$(1)_summary/selected_penalty_runs.gsea_jaccard.$(1).tsv: runs/$(1)_summary/summary.gsea_jaccard.$(1).tsv summarize_penalty_selection.py
	$(PYTHON) summarize_penalty_selection.py $$< runs/$(1)_summary/selected_penalty_runs.gsea_jaccard.$(1).tsv runs/$(1)_summary/penalty_selection.gsea_jaccard.$(1).pdf \
		--worst-tsv runs/$(1)_summary/worst_penalty_runs.gsea_jaccard.$(1).tsv \
		--clone $(1) \
		--overlap-col gsea_uniqueness \
		--overlap-label "GSEA uniqueness (1 - mean Jaccard)" \
		--overlap-higher-is-better

runs/$(1)_summary/worst_penalty_runs.gsea_jaccard.$(1).tsv runs/$(1)_summary/penalty_selection.gsea_jaccard.$(1).pdf: runs/$(1)_summary/selected_penalty_runs.gsea_jaccard.$(1).tsv
	@test -e $$@ || { echo "ERROR: expected output missing: $$@" >&2; exit 1; }

runs/$(1)_summary/selected_penalty_runs.top_genes_jaccard.$(1).tsv: runs/$(1)_summary/summary.top_genes_jaccard.$(1).tsv summarize_penalty_selection.py
	$(PYTHON) summarize_penalty_selection.py $$< runs/$(1)_summary/selected_penalty_runs.top_genes_jaccard.$(1).tsv runs/$(1)_summary/penalty_selection.top_genes_jaccard.$(1).pdf \
		--worst-tsv runs/$(1)_summary/worst_penalty_runs.top_genes_jaccard.$(1).tsv \
		--clone $(1) \
		--overlap-col top_gene_mean_pairwise_jaccard \
		--overlap-label "mean pairwise Jaccard of top up/down genes"

runs/$(1)_summary/worst_penalty_runs.top_genes_jaccard.$(1).tsv runs/$(1)_summary/penalty_selection.top_genes_jaccard.$(1).pdf: runs/$(1)_summary/selected_penalty_runs.top_genes_jaccard.$(1).tsv
	@test -e $$@ || { echo "ERROR: expected output missing: $$@" >&2; exit 1; }
endef

$(foreach clone,$(CLONES),$(eval $(call SUMMARY_TEMPLATE,$(clone))))

clean-metrics:
	rm -f $(L_CORR_SUMMARY_TSVS) $(GSEA_SUMMARY_TSVS) $(TOP_GENE_SUMMARY_TSVS) $(RUNTIME_SUMMARY_TSVS) $(RUN_METRICS) $(OLD_SUMMARY_OUTPUTS) $(OLD_CLONE_OUTPUTS)

clean-plots:
	rm -f $(L_CORR_PLOT_PDFS) $(GSEA_PLOT_PDFS) $(TOP_GENE_PLOT_PDFS) $(RUNTIME_PDFS) $(PENALTY_SELECTION_TSVS) $(PENALTY_WORST_TSVS) $(PENALTY_SELECTION_PDFS) $(PENALTY_GSEA_SELECTION_TSVS) $(PENALTY_GSEA_WORST_TSVS) $(PENALTY_GSEA_SELECTION_PDFS) $(PENALTY_TOP_GENE_SELECTION_TSVS) $(PENALTY_TOP_GENE_WORST_TSVS) $(PENALTY_TOP_GENE_SELECTION_PDFS) $(RUN_LOG_PDFS) $(RUN_FACTOR_HEATMAPS) $(RUN_FLOWTREES) $(RUN_FLOW_DISTANCES) $(OLD_PLOT_OUTPUTS) $(OLD_CLONE_OUTPUTS)

clean:
	rm -rf runs $(L_CORR_SUMMARY_TSVS) $(GSEA_SUMMARY_TSVS) $(TOP_GENE_SUMMARY_TSVS) $(RUNTIME_SUMMARY_TSVS) $(L_CORR_PLOT_PDFS) $(GSEA_PLOT_PDFS) $(TOP_GENE_PLOT_PDFS) $(RUNTIME_PDFS) $(PENALTY_SELECTION_TSVS) $(PENALTY_WORST_TSVS) $(PENALTY_SELECTION_PDFS) $(PENALTY_GSEA_SELECTION_TSVS) $(PENALTY_GSEA_WORST_TSVS) $(PENALTY_GSEA_SELECTION_PDFS) $(PENALTY_TOP_GENE_SELECTION_TSVS) $(PENALTY_TOP_GENE_WORST_TSVS) $(PENALTY_TOP_GENE_SELECTION_PDFS) $(RUN_LOG_PDFS) $(RUN_FACTOR_HEATMAPS) $(RUN_FLOWTREES) $(RUN_FLOW_DISTANCES) $(OLD_SUMMARY_OUTPUTS) $(OLD_PLOT_OUTPUTS) $(OLD_CLONE_OUTPUTS)
