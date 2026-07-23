export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
GEXFILTER := $(GEX_LINEAGE_DIR)/bin/gexFilter
GEXFACTOR := $(GEX_LINEAGE_DIR)/bin/gexFactor
BOLONDI_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
QUINN_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
PYTHON ?= /home/staklins/miniforge3/envs/pyplotting/bin/python

# CLONES can include TLS clones, such as tls1/tls2, and Quinn shorthand clones,
# such as quinn30 for Quinn CP 30.
CLONES ?= tls1 tls2 quinn4 quinn3 quinn2 quinn1

# Ability to downsample the input data for testing. Leave these empty to use the full data.
DOWNSAMPLE_CELLS ?= 
DOWNSAMPLE_GENES ?= 
DOWNSAMPLE_ENABLED := $(strip $(DOWNSAMPLE_CELLS)$(DOWNSAMPLE_GENES))
DOWNSAMPLE_CELLS_LABEL := $(if $(strip $(DOWNSAMPLE_CELLS)),$(DOWNSAMPLE_CELLS),all)
DOWNSAMPLE_GENES_LABEL := $(if $(strip $(DOWNSAMPLE_GENES)),$(DOWNSAMPLE_GENES),all)

# Sweep axes. The default grid varies only L1 and L-loading overlap.
# Set RIBO_MITO_VALUES="0 1" when you want that preprocessing toggle
# included in the same interpretable grid.
L_L1_STRENGTHS ?= 0 0.01 0.1 1 10
L_LOADING_OVERLAP_STRENGTHS ?= 0 0.01 0.1 1 10
RIBO_MITO_VALUES ?= 0

TREE_MODES ?= one
TREE_INDEX ?= 1
DIM ?= 5
DIM_VALUES ?= $(DIM)
SEED ?= 1
FILTER_TEST ?= moran

BASE_FILTER_OPTS := \
	--filter-test $(FILTER_TEST) \
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
$(if $(filter one,$(2)),runs/tree_inputs/$(1).one$(TREE_INDEX).nex,$(call TREE_INPUT,$(1)))
endef

define RIBO_MITO_OPTS
$(if $(filter 1 true yes,$(1)),--remove-ribo-mito-genes,)
endef

define GEX_INPUT
$(if $(DOWNSAMPLE_ENABLED),runs/gex_inputs/$(1).cells$(DOWNSAMPLE_CELLS_LABEL).genes$(DOWNSAMPLE_GENES_LABEL).tsv,$(call RAW_GEX_INPUT,$(1)))
endef

define CONDITION_NAME
T$(1)_Rm$(2)_K$(3)_Ll$(4)_Lo$(5)
endef

define CONDITION_OPTS
--L-l1-strength $(1) \
$(if $(filter-out 0 0.0 0.00 0.000,$(2)),--L-loading-overlap-strength $(2),)
endef

# Filtering is shared across factor-model parameter combinations. Include every
# setting that can alter the filtered matrix in its prefix.
define FILTER_NAME
T$(call TREE_MODE_NAME,$(2))_Rm$(3)_Ft$(subst -,_,$(FILTER_TEST))_S$(SEED)_C$(DOWNSAMPLE_CELLS_LABEL)_G$(DOWNSAMPLE_GENES_LABEL)
endef

define FILTER_PREFIX
runs/filtered/$(1).$(call FILTER_NAME,$(1),$(2),$(3))
endef

CONDITIONS := $(foreach tm,$(TREE_MODES),$(foreach rm,$(RIBO_MITO_VALUES),$(foreach k,$(DIM_VALUES),$(foreach ll1,$(L_L1_STRENGTHS),$(foreach lo,$(L_LOADING_OVERLAP_STRENGTHS),$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(tm)),$(rm),$(k),$(ll1),$(lo)))))))

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
SELECTION_TOP_N ?= 10
SELECTION_SUMMARY_TSVS := $(foreach clone,$(CLONES),runs/$(clone)_summary/selection_metrics.$(clone).tsv)
PENALTY_HEATMAP_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/penalty_heatmap.$(clone).pdf)
GSEA_PADJ_CUTOFF ?= 1.0
GSEA_TOP_N ?= 20
GSEA_HEATMAP_PADJ_CUTOFF ?= 0.05
GSEA_HEATMAP_PDFS := $(foreach clone,$(CLONES),runs/$(clone)_summary/gsea_heatmap.$(clone).pdf)
PAPER_SUMMARY_DIR := runs/summary
PAPER_TRADEOFF_PDF := $(PAPER_SUMMARY_DIR)/penalty_tradeoff.pdf
PAPER_SUMMARY_PDFS := $(PAPER_TRADEOFF_PDF)

all: run metrics plots
run: $(RUN_SUMMARIES)
metrics: $(SELECTION_SUMMARY_TSVS)
plots: $(RUN_LOG_PDFS) $(RUN_FACTOR_HEATMAPS) $(RUN_FLOWTREES) $(RUN_GSEAS) $(PENALTY_HEATMAP_PDFS) $(GSEA_HEATMAP_PDFS) $(PAPER_SUMMARY_PDFS)
paper-summary: $(PAPER_SUMMARY_PDFS)

print-config:
	@echo "CLONES=$(CLONES)"
	@echo "CONDITIONS=$(CONDITIONS)"
	@echo "L_L1_STRENGTHS=$(L_L1_STRENGTHS)"
	@echo "L_LOADING_OVERLAP_STRENGTHS=$(L_LOADING_OVERLAP_STRENGTHS)"
	@echo "TREE_MODES=$(TREE_MODES)"
	@echo "TREE_INDEX=$(TREE_INDEX)"
	@echo "RIBO_MITO_VALUES=$(RIBO_MITO_VALUES)"
	@echo "DIM_VALUES=$(DIM_VALUES)"
	@echo "DOWNSAMPLE_CELLS=$(DOWNSAMPLE_CELLS)"
	@echo "DOWNSAMPLE_GENES=$(DOWNSAMPLE_GENES)"
	@echo "FILTER_TEST=$(FILTER_TEST)"
	@echo "SELECTION_TOP_N=$(SELECTION_TOP_N)"
	@echo "GSEA_PADJ_CUTOFF=$(GSEA_PADJ_CUTOFF)"
	@echo "GSEA_TOP_N=$(GSEA_TOP_N)"
	@echo "GEXFILTER=$(GEXFILTER)"
	@echo "GEXFACTOR=$(GEXFACTOR)"
	@echo "BOLONDI_DIR=$(BOLONDI_DIR)"
	@echo "QUINN_DIR=$(QUINN_DIR)"
	@echo "PYPLOTTING_SIF=$(PYPLOTTING_SIF)"
	@echo "BASE_FILTER_OPTS=$(BASE_FILTER_OPTS)"

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

define ONE_TREE_INPUT_TEMPLATE
runs/tree_inputs/$(1).one$(TREE_INDEX).nex: $(call TREE_INPUT,$(1))
	mkdir -p runs/tree_inputs
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
		}' $$< > $$@
endef

$(if $(filter one,$(TREE_MODES)),$(foreach clone,$(CLONES),$(eval $(call ONE_TREE_INPUT_TEMPLATE,$(clone)))))

define FILTER_TEMPLATE
$(call FILTER_PREFIX,$(1),$(2),$(3)).filtered.tsv: \
$(call RUN_TREE_INPUT,$(1),$(2)) $(call GEX_INPUT,$(1))
	mkdir -p runs/filtered
	/usr/bin/time -o $(call FILTER_PREFIX,$(1),$(2),$(3)).time $(GEXFILTER) \
		--trees $(call RUN_TREE_INPUT,$(1),$(2)) \
		--expr $(call GEX_INPUT,$(1)) \
		--outprefix $(call FILTER_PREFIX,$(1),$(2),$(3)) \
		$(BASE_FILTER_OPTS) \
		$(call RIBO_MITO_OPTS,$(3)) > $(call FILTER_PREFIX,$(1),$(2),$(3)).term
endef

$(foreach clone,$(CLONES),$(foreach tm,$(TREE_MODES),$(foreach rm,$(RIBO_MITO_VALUES),$(eval $(call FILTER_TEMPLATE,$(clone),$(tm),$(rm))))))

define RUN_TEMPLATE
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).summary.tsv: \
$(call RUN_TREE_INPUT,$(1),$(2)) $(call FILTER_PREFIX,$(1),$(2),$(3)).filtered.tsv
	mkdir -p runs/$(1)
	/usr/bin/time -o runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).time $(GEXFACTOR) \
		--trees $(call RUN_TREE_INPUT,$(1),$(2)) \
		--expr $(call FILTER_PREFIX,$(1),$(2),$(3)).filtered.tsv \
		--outprefix runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)) \
		--dim $(4) \
		--seed $(SEED) \
		$(call CONDITION_OPTS,$(5),$(6)) > runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).term

runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).term \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).time \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).F.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).L.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).L.top_genes.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).X.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).log \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).latent_flow.tsv \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).L.pearson_correlation.tsv: \
runs/$(1)/$(call CONDITION_NAME,$(call TREE_MODE_NAME,$(2)),$(3),$(4),$(5),$(6)).summary.tsv
	@test -e $$@ || { echo "ERROR: expected output missing after gexFactor run: $$@" >&2; exit 1; }
endef

$(foreach clone,$(CLONES),$(foreach tm,$(TREE_MODES),$(foreach rm,$(RIBO_MITO_VALUES),$(foreach k,$(DIM_VALUES),$(foreach ll1,$(L_L1_STRENGTHS),$(foreach lo,$(L_LOADING_OVERLAP_STRENGTHS),$(eval $(call RUN_TEMPLATE,$(clone),$(tm),$(rm),$(k),$(ll1),$(lo)))))))))

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
runs/$(1)_summary/selection_metrics.$(1).tsv: $(foreach cond,$(CONDITIONS),runs/$(1)/$(cond).log runs/$(1)/$(cond).L.pearson_correlation.tsv runs/$(1)/$(cond).L.top_genes.tsv) summarize_selection_metrics.py
	mkdir -p runs/$(1)_summary
	$(PYTHON) summarize_selection_metrics.py --runs-dir runs --output $$@ \
		--clone $(1) \
		--top-n $(SELECTION_TOP_N) \
		--conditions $(CONDITIONS)

runs/$(1)_summary/penalty_heatmap.$(1).pdf: runs/$(1)_summary/selection_metrics.$(1).tsv plot_penalty_heatmap.py
	$(PYTHON) plot_penalty_heatmap.py $$< $$@ --clone $(1) --no-annot

runs/$(1)_summary/gsea_heatmap.$(1).pdf: $(RUN_GSEAS) plot_gsea_heatmap.py
	mkdir -p runs/$(1)_summary
	$(PYTHON) plot_gsea_heatmap.py runs $$@ --clone $(1) --padj-cutoff $(GSEA_HEATMAP_PADJ_CUTOFF) --no-annot
endef

$(foreach clone,$(CLONES),$(eval $(call SUMMARY_TEMPLATE,$(clone))))

$(PAPER_TRADEOFF_PDF): $(SELECTION_SUMMARY_TSVS) $(RUN_GSEAS) plot_penalty_summary.py
	mkdir -p $(PAPER_SUMMARY_DIR)
	$(PYTHON) plot_penalty_summary.py \
		--figure tradeoff \
		--output $@ \
		--runs-dir runs \
		--gsea-padj-cutoff $(GSEA_HEATMAP_PADJ_CUTOFF) \
		--selection-summaries $(SELECTION_SUMMARY_TSVS)

clean-metrics:
	rm -f $(SELECTION_SUMMARY_TSVS) $(RUN_METRICS) runs/*_summary/*.tsv

clean-plots:
	rm -f $(PENALTY_HEATMAP_PDFS) $(RUN_LOG_PDFS) $(RUN_FACTOR_HEATMAPS) $(RUN_FLOWTREES) $(RUN_FLOW_DISTANCES) $(PAPER_SUMMARY_PDFS) runs/*_summary/*.pdf

clean:
	rm -rf runs
