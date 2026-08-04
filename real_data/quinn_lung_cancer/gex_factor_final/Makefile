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
CURR_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_factor_final
PREPROCESSED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/preprocessed_data
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_filter
PCA_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/gex_pca
PCA_UMAP_DIR := $(BENCHMARKS_DIR)/real_data/quinn_lung_cancer/pca_umap
CLONES := 1 2 3 4

# Selected from ../gex_factor_dim_sweep/quinn_*.dim_summary.pdf.
DIM_1 := 1
DIM_2 := 1
DIM_3 := 1
DIM_4 := 1
NO_RIBO_MITO_DIM_1 := 10
NO_RIBO_MITO_DIM_2 := 10
NO_RIBO_MITO_DIM_3 := 10
NO_RIBO_MITO_DIM_4 := 10
CURRENT_DIM = $(shell test -f .chosen_dim.$(1) && sed -n '1p' .chosen_dim.$(1))
DIM_CHANGED = $(if $(filter $(DIM_$(1)),$(call CURRENT_DIM,$(1))),,FORCE)
NO_RIBO_MITO_DIM_CHANGED = $(if $(filter $(NO_RIBO_MITO_DIM_$(1)),$(call CURRENT_DIM,$(1)_no_ribo_mito)),,FORCE)

FIT_SUMMARIES := $(CLONES:%=quinn_%.fit.summary.tsv)
LOG_PLOTS := $(CLONES:%=quinn_%.fit.log.pdf)
PANEL_PLOTS := $(CLONES:%=quinn_%.fit.panel.png)
BROWNIAN_VARIANCE_PLOTS := $(CLONES:%=quinn_%.fit.brownian_variance.pdf)
GSEA := $(CLONES:%=quinn_%.fit.L.gsea.tsv)
FACTOR_HEATMAPS := $(CLONES:%=quinn_%.fit.L.factor_heatmap.pdf)
ACTIVE_FACTOR_HEATMAPS := $(CLONES:%=quinn_%.fit.L.active_factor_heatmap.pdf)
FLOW_TREES := $(CLONES:%=quinn_%.fit.latent_flow.factor_trees.pdf)
PCA_TREE_DEPTH := $(CLONES:%=quinn_%.fit.latent_flow.pca.tree_depth.pdf)
PCA_LATENT_DISTANCE := $(CLONES:%=quinn_%.fit.latent_flow.pca.latent_distance_from_root.pdf)
PCA_PLOT_STAMPS := $(CLONES:%=quinn_%.fit.latent_flow.pca.plots.done)
WADDINGTON_PLOTS := $(CLONES:%=quinn_%.fit.latent_flow.waddington.pdf)
NO_RIBO_MITO_FIT_SUMMARIES := $(CLONES:%=quinn_%_no_ribo_mito.fit.summary.tsv)
NO_RIBO_MITO_LOG_PLOTS := $(CLONES:%=quinn_%_no_ribo_mito.fit.log.pdf)
NO_RIBO_MITO_PANEL_PLOTS := $(CLONES:%=quinn_%_no_ribo_mito.fit.panel.png)
NO_RIBO_MITO_BROWNIAN_VARIANCE_PLOTS := $(CLONES:%=quinn_%_no_ribo_mito.fit.brownian_variance.pdf)
NO_RIBO_MITO_GSEA := $(CLONES:%=quinn_%_no_ribo_mito.fit.L.gsea.tsv)
NO_RIBO_MITO_FACTOR_HEATMAPS := $(CLONES:%=quinn_%_no_ribo_mito.fit.L.factor_heatmap.pdf)
NO_RIBO_MITO_ACTIVE_FACTOR_HEATMAPS := $(CLONES:%=quinn_%_no_ribo_mito.fit.L.active_factor_heatmap.pdf)
NO_RIBO_MITO_FLOW_TREES := $(CLONES:%=quinn_%_no_ribo_mito.fit.latent_flow.factor_trees.pdf)
NO_RIBO_MITO_PCA_TREE_DEPTH := $(CLONES:%=quinn_%_no_ribo_mito.fit.latent_flow.pca.tree_depth.pdf)
NO_RIBO_MITO_PCA_LATENT_DISTANCE := $(CLONES:%=quinn_%_no_ribo_mito.fit.latent_flow.pca.latent_distance_from_root.pdf)
NO_RIBO_MITO_PCA_PLOT_STAMPS := $(CLONES:%=quinn_%_no_ribo_mito.fit.latent_flow.pca.plots.done)
NO_RIBO_MITO_WADDINGTON_PLOTS := $(CLONES:%=quinn_%_no_ribo_mito.fit.latent_flow.waddington.pdf)
ACTIVE_FACTOR_CORRELATIONS := quinn_1_quinn_2_no_ribo_mito.fit.L.active_factor_correlations.pdf

.PHONY: all fit plotting archive clean clean-gexlineage FORCE
all: fit plotting
fit: $(FIT_SUMMARIES) $(NO_RIBO_MITO_FIT_SUMMARIES)
plotting: $(LOG_PLOTS) $(PANEL_PLOTS) $(BROWNIAN_VARIANCE_PLOTS) $(GSEA) $(FACTOR_HEATMAPS) $(ACTIVE_FACTOR_HEATMAPS) \
	$(FLOW_TREES) $(PCA_PLOT_STAMPS) $(WADDINGTON_PLOTS) \
	$(NO_RIBO_MITO_LOG_PLOTS) $(NO_RIBO_MITO_PANEL_PLOTS) $(NO_RIBO_MITO_BROWNIAN_VARIANCE_PLOTS) $(NO_RIBO_MITO_GSEA) \
	$(NO_RIBO_MITO_FACTOR_HEATMAPS) $(NO_RIBO_MITO_ACTIVE_FACTOR_HEATMAPS) \
	$(NO_RIBO_MITO_FLOW_TREES) $(NO_RIBO_MITO_PCA_PLOT_STAMPS) \
	$(NO_RIBO_MITO_WADDINGTON_PLOTS) $(ACTIVE_FACTOR_CORRELATIONS)

archive:
	@archive_dir=archive_$$(date +%Y-%m-%d_%H-%M-%S); \
	shopt -s nullglob; \
	outputs=( quinn_*.fit.* .chosen_dim .chosen_dim.* .fit_config.* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No final-factor outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

quinn_%.fit.summary.tsv: $(FILTERED_DIR)/quinn_%.tree.nex $(FILTERED_DIR)/quinn_%.moran.filtered.tsv $(PCA_DIR)/quinn_%.pca.eigenvectors.tsv $$(call DIM_CHANGED,$$*)
	/usr/bin/time -o quinn_$*.fit.time $(GEXFACTOR) --trees $(FILTERED_DIR)/quinn_$*.tree.nex --expr $(FILTERED_DIR)/quinn_$*.moran.filtered.tsv --outprefix quinn_$*.fit --dim $(DIM_$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_$*.fit.term
	@printf '%s\n' '$(DIM_$*)' > .chosen_dim.$*

quinn_%_no_ribo_mito.fit.summary.tsv: $(FILTERED_DIR)/quinn_%.tree.nex $(FILTERED_DIR)/quinn_%_no_ribo_mito.moran.filtered.tsv $(PCA_DIR)/quinn_%_no_ribo_mito.pca.eigenvectors.tsv $$(call NO_RIBO_MITO_DIM_CHANGED,$$*)
	/usr/bin/time -o quinn_$*_no_ribo_mito.fit.time $(GEXFACTOR) --trees $(FILTERED_DIR)/quinn_$*.tree.nex --expr $(FILTERED_DIR)/quinn_$*_no_ribo_mito.moran.filtered.tsv --outprefix quinn_$*_no_ribo_mito.fit --dim $(NO_RIBO_MITO_DIM_$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > quinn_$*_no_ribo_mito.fit.term
	@printf '%s\n' '$(NO_RIBO_MITO_DIM_$*)' > .chosen_dim.$*_no_ribo_mito

quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv quinn_%.fit.log quinn_%.fit.latent_flow.tsv: quinn_%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

quinn_%.fit.log.pdf: quinn_%.fit.summary.tsv
	@test -e quinn_$*.fit.log || { echo "ERROR: expected fit log missing: quinn_$*.fit.log" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_optimization_log.py /mnt/cwd/quinn_$*.fit.log -1 0 /mnt/cwd/quinn_$*.fit.log

quinn_%.fit.panel.png: quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv quinn_%.fit.summary.tsv $(SRC)/plot_fit_panel.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_fit_panel.py /mnt/cwd/quinn_$*.fit.F.tsv /mnt/cwd/quinn_$*.fit.L.tsv /mnt/cwd/quinn_$*.fit.X.tsv /mnt/cwd/quinn_$*.fit.summary.tsv /mnt/cwd/$@

quinn_%.fit.brownian_variance.pdf: quinn_%.fit.summary.tsv $(SRC)/plot_brownian_variance.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_brownian_variance.py /mnt/cwd/$< /mnt/cwd/$@

quinn_%.fit.L.gsea.tsv: quinn_%.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< human /mnt/cwd/quinn_$*.fit.L.gsea > quinn_$*.fit.gsea.term

quinn_%.fit.L.factor_heatmap.pdf: quinn_%.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_heatmap.py /mnt/cwd/$< /mnt/cwd/$@

quinn_%.fit.L.active_factor_heatmap.pdf: quinn_%.fit.L.tsv quinn_%.fit.summary.tsv $(SRC)/plot_factor_heatmap.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_heatmap.py /mnt/cwd/quinn_$*.fit.L.tsv /mnt/cwd/$@ \
		--fit-summary /mnt/cwd/quinn_$*.fit.summary.tsv

quinn_1_quinn_2_no_ribo_mito.fit.L.active_factor_correlations.pdf: \
		quinn_1_no_ribo_mito.fit.L.tsv quinn_1_no_ribo_mito.fit.summary.tsv \
		quinn_2_no_ribo_mito.fit.L.tsv quinn_2_no_ribo_mito.fit.summary.tsv \
		quinn_4_no_ribo_mito.fit.L.tsv quinn_4_no_ribo_mito.fit.summary.tsv \
		$(SRC)/plot_cross_sample_factor_correlations.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_cross_sample_factor_correlations.py \
		--comparison CP1 /mnt/cwd/quinn_1_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_1_no_ribo_mito.fit.summary.tsv CP2 /mnt/cwd/quinn_2_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_2_no_ribo_mito.fit.summary.tsv \
		--comparison CP1 /mnt/cwd/quinn_1_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_1_no_ribo_mito.fit.summary.tsv CP4 /mnt/cwd/quinn_4_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_4_no_ribo_mito.fit.summary.tsv \
		--comparison CP2 /mnt/cwd/quinn_2_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_2_no_ribo_mito.fit.summary.tsv CP4 /mnt/cwd/quinn_4_no_ribo_mito.fit.L.tsv /mnt/cwd/quinn_4_no_ribo_mito.fit.summary.tsv \
		--out-prefix /mnt/cwd/quinn_1_quinn_2_no_ribo_mito.fit.L.active_factor_correlations

quinn_%.fit.latent_flow.factor_trees.pdf: quinn_%.fit.latent_flow.tsv quinn_%.fit.summary.tsv $(SRC)/plot_latent_flow_trees.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_trees.py /mnt/cwd/quinn_$*.fit.latent_flow.tsv /mnt/cwd/quinn_$*.fit.latent_flow --fit-summary /mnt/cwd/quinn_$*.fit.summary.tsv

quinn_%_no_ribo_mito.fit.latent_flow.factor_trees.pdf: quinn_%_no_ribo_mito.fit.latent_flow.tsv quinn_%_no_ribo_mito.fit.summary.tsv $(SRC)/plot_latent_flow_trees.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_trees.py /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow.tsv /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow --fit-summary /mnt/cwd/quinn_$*_no_ribo_mito.fit.summary.tsv

quinn_%.fit.latent_flow.pca.plots.done: quinn_%.fit.latent_flow.tsv $(SRC)/plot_latent_flow_pca.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_pca.py /mnt/cwd/quinn_$*.fit.latent_flow.tsv /mnt/cwd/quinn_$*.fit.latent_flow
	@test -s quinn_$*.fit.latent_flow.pca.tree_depth.pdf
	@test -s quinn_$*.fit.latent_flow.pca.latent_distance_from_root.pdf
	@touch $@

quinn_%_no_ribo_mito.fit.latent_flow.pca.plots.done: quinn_%_no_ribo_mito.fit.latent_flow.tsv $(SRC)/plot_latent_flow_pca.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_pca.py /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow.tsv /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow
	@test -s quinn_$*_no_ribo_mito.fit.latent_flow.pca.tree_depth.pdf
	@test -s quinn_$*_no_ribo_mito.fit.latent_flow.pca.latent_distance_from_root.pdf
	@touch $@

quinn_%.fit.latent_flow.waddington.pdf: quinn_%.fit.latent_flow.tsv quinn_%.fit.summary.tsv $(PREPROCESSED_DIR)/quinn_%_barcode_tissue.tsv $(SRC)/plot_latent_flow_waddington.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_waddington.py /mnt/cwd/quinn_$*.fit.latent_flow.tsv /mnt/cwd/quinn_$*.fit.latent_flow --fit-summary /mnt/cwd/quinn_$*.fit.summary.tsv 
#--barcode-discrete-label-tsv /mnt/preprocessed_data/quinn_$*_barcode_tissue.tsv

quinn_%_no_ribo_mito.fit.latent_flow.waddington.pdf: quinn_%_no_ribo_mito.fit.latent_flow.tsv quinn_%_no_ribo_mito.fit.summary.tsv $(PREPROCESSED_DIR)/quinn_%_barcode_tissue.tsv $(SRC)/plot_latent_flow_waddington.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd --bind $(PREPROCESSED_DIR):/mnt/preprocessed_data $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_waddington.py /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow.tsv /mnt/cwd/quinn_$*_no_ribo_mito.fit.latent_flow --fit-summary /mnt/cwd/quinn_$*_no_ribo_mito.fit.summary.tsv 
#--barcode-discrete-label-tsv /mnt/preprocessed_data/quinn_$*_barcode_tissue.tsv

clean-gexlineage:
	rm -f quinn_*.fit.*

clean: clean-gexlineage
	rm -f .chosen_dim .chosen_dim.*
