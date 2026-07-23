export SHELL := /usr/bin/bash

.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:

BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := $(CONTAINERS)/pyplotting/pyplotting.sif
SCANPY_SIF := $(CONTAINERS)/scanpy/scanpy.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
GEXLINEAGE := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexLineage
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/final
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
CLONES := 1 2 3 4

# Choose this after inspecting ../dim_sweep/quinn_*.dim_summary.pdf.
DIM ?=
ifeq ($(strip $(DIM)),)
$(error DIM is required; run make DIM=<chosen_dimension>)
endif

CURRENT_DIM := $(shell test -f .chosen_dim && sed -n '1p' .chosen_dim)
ifneq ($(strip $(DIM)),$(strip $(CURRENT_DIM)))
DIM_CHANGED := FORCE
endif

FIT_SUMMARIES := $(CLONES:%=quinn_%.fit.summary.tsv)
ORIG_EMBEDDINGS := $(CLONES:%=quinn_%_gex.umap.pdf)
LOG_PLOTS := $(CLONES:%=quinn_%.fit.log.pdf)
PANEL_PLOTS := $(CLONES:%=quinn_%.fit.panel.png)
GSEA := $(CLONES:%=quinn_%.fit.L.gsea.tsv)
FACTOR_HEATMAPS := $(CLONES:%=quinn_%.fit.L.factor_heatmap.pdf)
FLOW_TREES := $(CLONES:%=quinn_%.fit.latent_flow.factor_trees.pdf)
EUCLIDEAN_PSEUDOTIME := $(CLONES:%=quinn_%.fit.latent_flow.euclidean_factor_pseudotime.pdf)

.PHONY: all fit plotting clean-gexlineage FORCE
all: fit plotting
fit: $(FIT_SUMMARIES)
plotting: $(ORIG_EMBEDDINGS) $(LOG_PLOTS) $(PANEL_PLOTS) $(GSEA) $(FACTOR_HEATMAPS) $(FLOW_TREES) $(EUCLIDEAN_PSEUDOTIME)

quinn_%.fit.summary.tsv: $(PREPROCESSED_DIR)/quinn_%.vine.nex $(PREPROCESSED_DIR)/quinn_%_gex.tsv $(DIM_CHANGED)
	NTAXA=$$(wc -l < $(PREPROCESSED_DIR)/quinn_$*_gex.tsv); if [ $$NTAXA -gt 1000 ]; then FILTERTEST=moran; else FILTERTEST=lrt; fi; \
	/usr/bin/time -o quinn_$*.fit.time $(GEXLINEAGE) --trees $(PREPROCESSED_DIR)/quinn_$*.vine.nex --expr $(PREPROCESSED_DIR)/quinn_$*_gex.tsv --tree-total-time 1.0 --outprefix quinn_$*.fit --filter-test $$FILTERTEST --dim $(DIM) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --tree-index 1 --nthreads 1 --seed 1 > quinn_$*.fit.term
	@printf '%s\n' '$(DIM)' > .chosen_dim

quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv quinn_%.fit.log quinn_%.fit.latent_flow.tsv: quinn_%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexLineage run: $@" >&2; exit 1; }

quinn_%_gex.pca.pdf quinn_%_gex.umap.pdf quinn_%_gex.pca.tsv quinn_%_gex.umap.tsv: $(PREPROCESSED_DIR)/quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(SCANPY_SIF) python /mnt/src/scanpy_dim_reduction.py /mnt/preprocessed_data/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

quinn_%.fit.log.pdf: quinn_%.fit.summary.tsv
	@test -e quinn_$*.fit.log || { echo "ERROR: expected fit log missing: quinn_$*.fit.log" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_optimization_log.py /mnt/cwd/quinn_$*.fit.log -1 0 /mnt/cwd/quinn_$*.fit.log

quinn_%.fit.panel.png: quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_fit_panel.py /mnt/cwd/quinn_$*.fit.F.tsv /mnt/cwd/quinn_$*.fit.L.tsv /mnt/cwd/quinn_$*.fit.X.tsv /mnt/cwd/$@

quinn_%.fit.L.gsea.tsv: quinn_%.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< human /mnt/cwd/quinn_$*.fit.L.gsea > quinn_$*.fit.gsea.term

quinn_%.fit.L.factor_heatmap.pdf: quinn_%.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_heatmap.py /mnt/cwd/$< /mnt/cwd/$@

quinn_%.fit.latent_flow.factor_trees.pdf quinn_%.fit.latent_flow.total_euclidean_distance.tsv: quinn_%.fit.summary.tsv $(PREPROCESSED_DIR)/quinn_%_barcode_tissue.tsv
	@test -e quinn_$*.fit.latent_flow.tsv || { echo "ERROR: expected latent-flow file missing: quinn_$*.fit.latent_flow.tsv" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow.py /mnt/cwd/quinn_$*.fit.latent_flow.tsv /mnt/cwd/quinn_$*.fit.latent_flow --barcode-discrete-label-tsv /mnt/preprocessed_data/quinn_$*_barcode_tissue.tsv

quinn_%.fit.latent_flow.euclidean_factor_pseudotime.pdf: quinn_%.fit.latent_flow.total_euclidean_distance.tsv quinn_%_gex.pca.tsv quinn_%_gex.umap.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_euclidean_distance_pseudotime.py /mnt/cwd/quinn_$*_gex.pca.tsv /mnt/cwd/quinn_$*_gex.umap.tsv /mnt/cwd/quinn_$*.fit.latent_flow.total_euclidean_distance.tsv /mnt/cwd/quinn_$*.fit.latent_flow

clean-gexlineage:
	rm -f quinn_*.fit.*
