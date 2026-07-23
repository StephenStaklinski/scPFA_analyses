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
GEXFACTOR := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexFactor

CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/factor_final
SHARED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/preprocessed_data
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_filter
SAMPLES := tls1 tls2

# Set on the command line after inspecting ../factor_dim_sweep/*.dim_summary.pdf,
# for example: make DIM=8
DIM ?=
ifeq ($(strip $(DIM)),)
$(error DIM is required; run make DIM=<chosen_dimension>)
endif

CURRENT_DIM := $(shell test -f .chosen_dim && sed -n '1p' .chosen_dim)
ifneq ($(strip $(DIM)),$(strip $(CURRENT_DIM)))
DIM_CHANGED := FORCE
endif

FIT_SUMMARIES := $(SAMPLES:%=%.fit.summary.tsv)
ORIG_EMBEDDINGS := $(SAMPLES:%=%_gex.umap.pdf)
LOG_PLOTS := $(SAMPLES:%=%.fit.log.pdf)
PANEL_PLOTS := $(SAMPLES:%=%.fit.panel.png)
GSEA := $(SAMPLES:%=%.fit.L.gsea.tsv)
FACTOR_HEATMAPS := $(SAMPLES:%=%.fit.L.factor_heatmap.pdf)
FLOW_TREES := $(SAMPLES:%=%.fit.latent_flow.factor_trees.pdf)
EUCLIDEAN_PSEUDOTIME := $(SAMPLES:%=%.fit.latent_flow.euclidean_factor_pseudotime.pdf)

.PHONY: all fit plotting clean-gexlineage FORCE
all: fit plotting
fit: $(FIT_SUMMARIES)
plotting: $(ORIG_EMBEDDINGS) $(LOG_PLOTS) $(PANEL_PLOTS) $(GSEA) $(FACTOR_HEATMAPS) $(FLOW_TREES) $(EUCLIDEAN_PSEUDOTIME)

%.fit.summary.tsv: $(FILTERED_DIR)/%.tree.nex $(FILTERED_DIR)/%.moran.filtered.tsv $(DIM_CHANGED)
	/usr/bin/time -o $*.fit.time $(GEXFACTOR) --trees $(FILTERED_DIR)/$*.tree.nex --expr $(FILTERED_DIR)/$*.moran.filtered.tsv --outprefix $*.fit --dim $(DIM) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > $*.fit.term
	@printf '%s\n' '$(DIM)' > .chosen_dim

%.fit.F.tsv %.fit.L.tsv %.fit.X.tsv %.fit.log %.fit.latent_flow.tsv: %.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

%_gex.pca.pdf %_gex.umap.pdf %_gex.pca.tsv %_gex.umap.tsv: $(SHARED_DIR)/%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(SHARED_DIR):/mnt/shared $(SCANPY_SIF) python /mnt/src/scanpy_dim_reduction.py /mnt/shared/$*_gex.tsv /mnt/cwd/$*_gex

%.fit.log.pdf: %.fit.summary.tsv
	@test -e $*.fit.log || { echo "ERROR: expected fit log missing: $*.fit.log" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_optimization_log.py /mnt/cwd/$*.fit.log -1 0 /mnt/cwd/$*.fit.log

%.fit.panel.png: %.fit.F.tsv %.fit.L.tsv %.fit.X.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_fit_panel.py /mnt/cwd/$*.fit.F.tsv /mnt/cwd/$*.fit.L.tsv /mnt/cwd/$*.fit.X.tsv /mnt/cwd/$@

%.fit.L.gsea.tsv: %.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< mouse /mnt/cwd/$*.fit.L.gsea > $*.fit.gsea.term

%.fit.L.factor_heatmap.pdf: %.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_heatmap.py /mnt/cwd/$< /mnt/cwd/$@

%.fit.latent_flow.factor_trees.pdf %.fit.latent_flow.total_euclidean_distance.tsv: %.fit.summary.tsv
	@test -e $*.fit.latent_flow.tsv || { echo "ERROR: expected latent-flow file missing: $*.fit.latent_flow.tsv" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow.py /mnt/cwd/$*.fit.latent_flow.tsv /mnt/cwd/$*.fit.latent_flow

%.fit.latent_flow.euclidean_factor_pseudotime.pdf: %.fit.latent_flow.total_euclidean_distance.tsv %_gex.pca.tsv %_gex.umap.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_euclidean_distance_pseudotime.py /mnt/cwd/$*_gex.pca.tsv /mnt/cwd/$*_gex.umap.tsv /mnt/cwd/$*.fit.latent_flow.total_euclidean_distance.tsv /mnt/cwd/$*.fit.latent_flow

clean-gexlineage:
	rm -f *.fit.*
