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

CURR_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_factor_final
FILTERED_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_filter
PCA_DIR := $(BENCHMARKS_DIR)/real_data/bolondi_tls_embryogenesis/gex_pca

SAMPLES := tls1 tls2

# Selected from ../gex_factor_dim_sweep/*.dim_summary.pdf.
DIM_tls1 := 10
DIM_tls2 := 10
NO_RIBO_MITO_DIM_tls1 := 10
NO_RIBO_MITO_DIM_tls2 := 10
FIT_CONFIG = dim=$(DIM_$(1)) latent_variance=fixed_one
NO_RIBO_MITO_FIT_CONFIG = dim=$(NO_RIBO_MITO_DIM_$(1)) latent_variance=fixed_one
CURRENT_FIT_CONFIG = $(shell test -f .fit_config.$(1) && sed -n '1p' .fit_config.$(1))
FIT_CONFIG_CHANGED = $(if $(shell test '$(call CURRENT_FIT_CONFIG,$(1))' = '$(call FIT_CONFIG,$(1))' && printf yes),,FORCE)
NO_RIBO_MITO_FIT_CONFIG_CHANGED = $(if $(shell test '$(call CURRENT_FIT_CONFIG,$(1)_no_ribo_mito)' = '$(call NO_RIBO_MITO_FIT_CONFIG,$(1))' && printf yes),,FORCE)

FIT_SUMMARIES := $(SAMPLES:%=%.fit.summary.tsv)
LOG_PLOTS := $(SAMPLES:%=%.fit.log.pdf)
PANEL_PLOTS := $(SAMPLES:%=%.fit.panel.png)
BROWNIAN_VARIANCE_PLOTS := $(SAMPLES:%=%.fit.brownian_variance.pdf)
GSEA := $(SAMPLES:%=%.fit.L.gsea.tsv)
FACTOR_HEATMAPS := $(SAMPLES:%=%.fit.L.factor_heatmap.pdf)
ACTIVE_FACTOR_HEATMAPS := $(SAMPLES:%=%.fit.L.active_factor_heatmap.pdf)
FLOW_TREES := $(SAMPLES:%=%.fit.latent_flow.factor_trees.pdf)
WADDINGTON_PLOTS := $(SAMPLES:%=%.fit.latent_flow.waddington.pdf)
NO_RIBO_MITO_SAMPLES := $(SAMPLES:%=%_no_ribo_mito)
NO_RIBO_MITO_FIT_SUMMARIES := $(NO_RIBO_MITO_SAMPLES:%=%.fit.summary.tsv)
NO_RIBO_MITO_LOG_PLOTS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.log.pdf)
NO_RIBO_MITO_PANEL_PLOTS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.panel.png)
NO_RIBO_MITO_BROWNIAN_VARIANCE_PLOTS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.brownian_variance.pdf)
NO_RIBO_MITO_GSEA := $(NO_RIBO_MITO_SAMPLES:%=%.fit.L.gsea.tsv)
NO_RIBO_MITO_FACTOR_HEATMAPS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.L.factor_heatmap.pdf)
NO_RIBO_MITO_ACTIVE_FACTOR_HEATMAPS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.L.active_factor_heatmap.pdf)
NO_RIBO_MITO_FLOW_TREES := $(NO_RIBO_MITO_SAMPLES:%=%.fit.latent_flow.factor_trees.pdf)
NO_RIBO_MITO_WADDINGTON_PLOTS := $(NO_RIBO_MITO_SAMPLES:%=%.fit.latent_flow.waddington.pdf)
ACTIVE_FACTOR_CORRELATIONS := tls1_tls2_no_ribo_mito.fit.L.active_factor_correlations.pdf

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
	outputs=( *.fit.* .fit_config.* .chosen_dim .chosen_dim.* ); \
	if (( $${#outputs[@]} == 0 )); then \
		echo "No final-factor outputs to archive."; \
	else \
		mkdir -p "$$archive_dir"; \
		mv "$${outputs[@]}" "$$archive_dir"/; \
		echo "Archived $${#outputs[@]} outputs in $$archive_dir"; \
	fi

# $(PCA_DIR)/%.pca.eigenvectors.tsv
%.fit.summary.tsv: $(FILTERED_DIR)/%.tree.nex $(FILTERED_DIR)/%.moran.filtered.tsv $$(call FIT_CONFIG_CHANGED,$$*)
	/usr/bin/time -o $*.fit.time $(GEXFACTOR) --trees $(FILTERED_DIR)/$*.tree.nex --expr $(FILTERED_DIR)/$*.moran.filtered.tsv --outprefix $*.fit --dim $(DIM_$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > $*.fit.term
	@printf '%s\n' '$(call FIT_CONFIG,$*)' > .fit_config.$*

# $(PCA_DIR)/%_no_ribo_mito.pca.eigenvectors.tsv
%_no_ribo_mito.fit.summary.tsv: $(FILTERED_DIR)/%.tree.nex $(FILTERED_DIR)/%_no_ribo_mito.moran.filtered.tsv $$(call NO_RIBO_MITO_FIT_CONFIG_CHANGED,$$*)
	/usr/bin/time -o $*_no_ribo_mito.fit.time $(GEXFACTOR) --trees $(FILTERED_DIR)/$*.tree.nex --expr $(FILTERED_DIR)/$*_no_ribo_mito.moran.filtered.tsv --outprefix $*_no_ribo_mito.fit --dim $(NO_RIBO_MITO_DIM_$*) --L-l1-strength 0.1 --L-loading-overlap-strength 1 --seed 1 > $*_no_ribo_mito.fit.term
	@printf '%s\n' '$(call NO_RIBO_MITO_FIT_CONFIG,$*)' > .fit_config.$*_no_ribo_mito

%.fit.F.tsv %.fit.L.tsv %.fit.X.tsv %.fit.log %.fit.latent_flow.tsv: | %.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexFactor run: $@" >&2; exit 1; }

%.fit.log.pdf: %.fit.summary.tsv
	@test -e $*.fit.log || { echo "ERROR: expected fit log missing: $*.fit.log" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_optimization_log.py /mnt/cwd/$*.fit.log -1 0 /mnt/cwd/$*.fit.log

%.fit.panel.png: %.fit.F.tsv %.fit.L.tsv %.fit.X.tsv %.fit.summary.tsv $(SRC)/plot_fit_panel.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_fit_panel.py /mnt/cwd/$*.fit.F.tsv /mnt/cwd/$*.fit.L.tsv /mnt/cwd/$*.fit.X.tsv /mnt/cwd/$*.fit.summary.tsv /mnt/cwd/$@

%.fit.brownian_variance.pdf: %.fit.summary.tsv $(SRC)/plot_brownian_variance.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_brownian_variance.py /mnt/cwd/$< /mnt/cwd/$@

%.fit.L.gsea.tsv: %.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) /mnt/src/run_factor_gsea.R /mnt/cwd/$< mouse /mnt/cwd/$*.fit.L.gsea > $*.fit.gsea.term

%.fit.L.factor_heatmap.pdf: %.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_factor_heatmap.py /mnt/cwd/$< /mnt/cwd/$@

%.fit.L.active_factor_heatmap.pdf: %.fit.L.tsv %.fit.summary.tsv $(SRC)/plot_factor_heatmap.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_heatmap.py /mnt/cwd/$*.fit.L.tsv /mnt/cwd/$@ \
		--fit-summary /mnt/cwd/$*.fit.summary.tsv

tls1_tls2_no_ribo_mito.fit.L.active_factor_correlations.pdf: \
		tls1_no_ribo_mito.fit.L.tsv tls1_no_ribo_mito.fit.summary.tsv \
		tls2_no_ribo_mito.fit.L.tsv tls2_no_ribo_mito.fit.summary.tsv \
		$(SRC)/plot_cross_sample_factor_correlations.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) \
		python /mnt/src/plot_cross_sample_factor_correlations.py \
		--row-sample TLS1 /mnt/cwd/tls1_no_ribo_mito.fit.L.tsv /mnt/cwd/tls1_no_ribo_mito.fit.summary.tsv \
		--column-sample TLS2 /mnt/cwd/tls2_no_ribo_mito.fit.L.tsv /mnt/cwd/tls2_no_ribo_mito.fit.summary.tsv \
		--out-prefix /mnt/cwd/tls1_tls2_no_ribo_mito.fit.L.active_factor_correlations

%.fit.latent_flow.factor_trees.pdf: %.fit.latent_flow.tsv %.fit.summary.tsv $(SRC)/plot_latent_flow_trees.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_trees.py /mnt/cwd/$*.fit.latent_flow.tsv /mnt/cwd/$*.fit.latent_flow --fit-summary /mnt/cwd/$*.fit.summary.tsv

%.fit.latent_flow.waddington.pdf: %.fit.latent_flow.tsv %.fit.summary.tsv $(SRC)/plot_latent_flow_waddington.py
	singularity exec --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd $(PYPLOTTING_SIF) python /mnt/src/plot_latent_flow_waddington.py /mnt/cwd/$*.fit.latent_flow.tsv /mnt/cwd/$*.fit.latent_flow --fit-summary /mnt/cwd/$*.fit.summary.tsv

clean-gexlineage:
	rm -f *.fit.*

clean: clean-gexlineage
	rm -f .chosen_dim .chosen_dim.* .fit_config.*
