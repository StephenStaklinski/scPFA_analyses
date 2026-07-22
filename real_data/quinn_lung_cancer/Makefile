export SHELL=/usr/bin/bash

QUINNROOT := ../data

# Prevent Make from removing files it thinks it should clean up
.NOTINTERMEDIATE:
.SECONDARY:
.PRECIOUS:
.SECONDEXPANSION:


BENCHMARKS_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks
SRC := $(BENCHMARKS_DIR)/src
CONTAINERS := $(BENCHMARKS_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
PYPLOTTING_SIF := ${CONTAINERS}/pyplotting/pyplotting.sif
SCANPY_SIF := ${CONTAINERS}/scanpy/scanpy.sif
CLUSTERPROFILER_SIF := $(CONTAINERS)/clusterprofiler/clusterprofiler.sif
MONOCLE_SIF := ${CONTAINERS}/monocle/monocle.sif
VINE := /home/staklins/projects/vine_project/vine/bin/vine
GEXLINEAGE := /home/staklins/projects/gex_lineage_project/gex_lineage/bin/gexLineage

CURR_DIR := $(BENCHMARKS_DIR)/quinn_lung_cancer

DATA_DIR_CHAR := ./data/GSE161363/GSM4905334_character_matrices_5k
DATA_DIR_GEX := ./data/GSE161363/GSM4905335_barcodes_genes_matrix_5k


# Excluded in the original paper barcode quality control: 5 16 18 25 33 38 39 41 50 53 65 69 75 81 87 88 93
# Ignore because they have <1000 cells: 6 7 8 9 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 28 29 30 31 32 34 35 36 37 40 42 43 44 45 46 47 48 49 51 52 54 55 56 57 58 59 60 61 62 63 64 66 67 68 70 71 72 73 74 76 77 78 79 80 82 83 84 85 86 89 90 91 92 94 95 96 97 98 99 100
# Failed (due to memory): 1
QUINNCPS := 1 2 3 4  

# Dimension sweep settings
DIM_VALUES ?= 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 30 40 50

FIT_PREFIXES := $(foreach cp,$(QUINNCPS),$(foreach dim,$(DIM_VALUES),quinn_$(cp).K$(dim).fit))

FIT_CP = $(firstword $(subst .K, ,$(1)))
FIT_DIM = $(word 2,$(subst .K, ,$(1)))
FIT_TREE = quinn_$(call FIT_CP,$(1)).vine.nex
FIT_GEX = quinn_$(call FIT_CP,$(1))_gex.tsv
FIT_BARCODE_TISSUE = quinn_$(call FIT_CP,$(1))_barcode_tissue.tsv
FIT_PCA = quinn_$(call FIT_CP,$(1))_gex.pca.tsv
FIT_UMAP = quinn_$(call FIT_CP,$(1))_gex.umap.tsv
DIM_FIT_LOGS = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.log)
DIM_FIT_SUMMARIES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.summary.tsv)
DIM_FIT_GSEAS = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.L.gsea.tsv)
DIM_FIT_TIMES = $(foreach dim,$(DIM_VALUES),quinn_$(1).K$(dim).fit.time)


CHARMATRICES := $(QUINNCPS:%=quinn_%_character_matrix.tsv)
BARCODE_TISSUES := $(QUINNCPS:%=quinn_%_barcode_tissue.tsv)
GEXMATS := $(QUINNCPS:%=quinn_%_gex.tsv)
VINETREES := $(QUINNCPS:%=quinn_%.vine.nex)
GEXLINEAGERESULTS := $(FIT_PREFIXES:%=%.summary.tsv)
LOGPLOTS := $(FIT_PREFIXES:%=%.log.pdf)
PANELPLOTS := $(FIT_PREFIXES:%=%.panel.png)
GEXLINEAGEGSEA := $(FIT_PREFIXES:%=%.L.gsea.tsv)
FACTORHEATMAPS := $(FIT_PREFIXES:%=%.L.factor_heatmap.pdf)
DIMSUMMARYPLOTS := $(QUINNCPS:%=quinn_%.dim_summary.pdf)
ORIGEMBEDDINGS := $(QUINNCPS:%=quinn_%_gex.umap.pdf)
PALANTIRTRAJECTORIES := $(QUINNCPS:%=quinn_%_gex.palantir.pseudotime.tsv)
MONOCLETRAJECTORIES := $(QUINNCPS:%=quinn_%_gex.monocle.pseudotime.pdf)
FLOWTREES := $(FIT_PREFIXES:%=%.latent_flow.factor_trees.pdf)
EUCLIDEANPSEUDOTIME := $(FIT_PREFIXES:%=%.latent_flow.euclidean_factor_pseudotime.pdf)

all: prep_data run_vine gexlineage plotting dim_plots

prep_data: $(CHARMATRICES) $(BARCODE_TISSUES) $(GEXMATS) $(ORIGEMBEDDINGS)
palantir: $(PALANTIRTRAJECTORIES)
monocle: $(MONOCLETRAJECTORIES)
run_vine: $(VINETREES)
gexlineage: $(GEXLINEAGERESULTS)
plotting: $(LOGPLOTS) $(PANELPLOTS) $(FACTORHEATMAPS) $(GEXLINEAGEGSEA) $(FLOWTREES) $(EUCLIDEANPSEUDOTIME)
dim_plots: $(DIMSUMMARYPLOTS)


# Copy and reformat the provided character matrices
quinn_%_character_matrix.tsv:
	awk 'BEGIN{FS=OFS="\t"} { for(i=2;i<=NF;i++) if($$i=="-") $$i="-1"; print }' $(DATA_DIR_CHAR)/m5k_lg$*_character_matrix.alleleThresh.txt > $@

# Extract the tissue prefix from cell names such as RW.AAAC... into clone-specific metadata.
quinn_%_barcode_tissue.tsv: quinn_%_character_matrix.tsv
	awk 'BEGIN{FS=OFS="\t"; print "barcode", "tissue"} FNR > 1 { tissue=$$1; sub(/\..*$$/, "", tissue); print $$1, tissue }' $< > $@

# Subset the gene expression data to the cells present in the character matrix
quinn_%_gex.tsv: quinn_%_character_matrix.tsv $(DATA_DIR_GEX)/GSM4905335_matrix.5k.mtx $(DATA_DIR_GEX)/GSM4905335_barcodes.5k.tsv $(DATA_DIR_GEX)/GSM4905335_genes.5k.tsv
		singularity exec --bind $(DATA_DIR_GEX):/mnt/gex --bind $(SRC):/mnt/src --bind $(CURR_DIR):/mnt/cwd/ $(CASSIOPEIA_SIF) \
			python /mnt/src/subset_mtx_to_char_matrix_cells.py \
			/mnt/cwd/quinn_$*_character_matrix.tsv /mnt/gex/GSM4905335_matrix.5k.mtx \
			/mnt/gex/GSM4905335_barcodes.5k.tsv /mnt/gex/GSM4905335_genes.5k.tsv \
			/mnt/cwd/quinn_$*_gex.tsv > quinn_$*_gex.term

# Plot the original gex embeddings
quinn_%_gex.pca.pdf quinn_%_gex.umap.pdf quinn_%_gex.pca.tsv quinn_%_gex.umap.tsv: quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(SCANPY_SIF) \
		python /mnt/src/scanpy_dim_reduction.py /mnt/cwd/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

# Get gex trajectories with palantir
quinn_%_gex.palantir.pseudotime.tsv: quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(SCANPY_SIF) \
		python /mnt/src/palantir_trajectory.py /mnt/cwd/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

# Get gex trajectories with monocle
quinn_%_gex.monocle.pseudotime.pdf: quinn_%_gex.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(MONOCLE_SIF) \
		Rscript /mnt/src/monocle_trajectory.R /mnt/cwd/quinn_$*_gex.tsv /mnt/cwd/quinn_$*_gex

# Run VINE
quinn_%.vine.log quinn_%.vine.nwk quinn_%.vine.time quinn_%.vine.mean.nwk: quinn_%_character_matrix.tsv
	/usr/bin/time -o quinn_$*.vine.time $(VINE) -M 400 -c 100 -v 0 \
	 -i CRISPR quinn_$*_character_matrix.tsv --logf quinn_$*.vine.log --mean quinn_$*.vine.mean.nwk > quinn_$*.vine.nwk

# Convert nwk to nexus
quinn_%.vine.nex: quinn_%.vine.nwk
	singularity exec --bind $(CURR_DIR):/mnt/cwd/ --bind $(SRC):/mnt/src $(CASSIOPEIA_SIF) \
		python /mnt/src/nwk2nex.py /mnt/cwd/quinn_$*.vine.nwk /mnt/cwd/quinn_$*.vine.nex

# Fit the gene expression model
quinn_%.fit.summary.tsv: $$(call FIT_TREE,$$*) $$(call FIT_GEX,$$*)
	NTAXA=$$(wc -l < $(call FIT_GEX,$*)); \
	if [ $$NTAXA -gt 1000 ]; then \
		FILTERTEST=moran; \
	else \
		FILTERTEST=lrt; \
	fi; \
	/usr/bin/time -o quinn_$*.fit.time ${GEXLINEAGE} \
	--trees $(call FIT_TREE,$*) \
	--expr $(call FIT_GEX,$*) \
	--tree-total-time 1.0 \
	--outprefix quinn_$*.fit \
	--filter-test $$FILTERTEST \
	--dim $(call FIT_DIM,$*) \
	--L-l1-strength 0.1 \
	--L-loading-overlap-strength 1 \
	--tree-index 1 \
	--nthreads 1 \
	--seed 1 > quinn_$*.fit.term

quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv quinn_%.fit.log quinn_%.fit.latent_flow.tsv quinn_%.fit.time: quinn_%.fit.summary.tsv
	@test -e $@ || { echo "ERROR: expected output missing after gexLineage run: $@" >&2; exit 1; }

# Plot the fit optimization log results.  The log is written before the
# summary during a fit, so use the summary as the freshness prerequisite.
quinn_%.fit.log.pdf: quinn_%.fit.summary.tsv
	@test -e quinn_$*.fit.log || { echo "ERROR: expected fit log missing: quinn_$*.fit.log" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_optimization_log.py /mnt/cwd/quinn_$*.fit.log -1 0 /mnt/cwd/quinn_$*.fit.log

# Plot dimensionality-sweep diagnostics (one combined PDF)
quinn_%.dim_summary.pdf: $$(call DIM_FIT_SUMMARIES,$$*) $$(call DIM_FIT_GSEAS,$$*) $$(call DIM_FIT_TIMES,$$*) $(SRC)/plot_dim_log_likelihood.py
	@for log in $(call DIM_FIT_LOGS,$*); do \
		test -e $$log || { echo "ERROR: expected fit log missing: $$log" >&2; exit 1; }; \
	done
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_dim_log_likelihood.py \
			/mnt/cwd/quinn_$*.dim_summary.pdf \
			$(foreach dim,$(DIM_VALUES),/mnt/cwd/quinn_$*.K$(dim).fit.log)

# Plot the fit results
quinn_%.fit.panel.png: quinn_%.fit.F.tsv quinn_%.fit.L.tsv quinn_%.fit.X.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_fit_panel.py \
			/mnt/cwd/quinn_$*.fit.F.tsv \
			/mnt/cwd/quinn_$*.fit.L.tsv \
			/mnt/cwd/quinn_$*.fit.X.tsv \
			/mnt/cwd/quinn_$*.fit.panel.png

# Run GO BP GSEA on fitted L matrix
quinn_%.fit.L.gsea.tsv: quinn_%.fit.L.tsv
	singularity run --bind $(CURR_DIR):/mnt/cwd --bind $(SRC):/mnt/src $(CLUSTERPROFILER_SIF) \
		/mnt/src/run_factor_gsea.R /mnt/cwd/quinn_$*.fit.L.tsv human /mnt/cwd/quinn_$*.fit.L.gsea > quinn_$*.fit.gsea.term

# Plot the factor heatmap
quinn_%.fit.L.factor_heatmap.pdf: quinn_%.fit.L.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_heatmap.py \
			/mnt/cwd/quinn_$*.fit.L.tsv \
			/mnt/cwd/quinn_$*.fit.L.factor_heatmap.pdf

# Plot the latent flow results.  As above, latent_flow.tsv is an earlier
# side output of the fit, so let the summary drive freshness.
quinn_%.fit.latent_flow.factor_trees.pdf quinn_%.fit.latent_flow.total_euclidean_distance.tsv: quinn_%.fit.summary.tsv $$(call FIT_BARCODE_TISSUE,$$*)
	@test -e quinn_$*.fit.latent_flow.tsv || { echo "ERROR: expected latent-flow file missing: quinn_$*.fit.latent_flow.tsv" >&2; exit 1; }
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_latent_flow.py \
			/mnt/cwd/quinn_$*.fit.latent_flow.tsv \
			/mnt/cwd/quinn_$*.fit.latent_flow \
			--barcode-discrete-label-tsv /mnt/cwd/$(call FIT_BARCODE_TISSUE,$*)

# Plot euclidean latent-factor pseudotime on gex-only PCA/UMAP embeddings
quinn_%.fit.latent_flow.euclidean_factor_pseudotime.pdf: quinn_%.fit.latent_flow.total_euclidean_distance.tsv $$(call FIT_PCA,$$*) $$(call FIT_UMAP,$$*)
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_factor_euclidean_distance_pseudotime.py \
			/mnt/cwd/$(call FIT_PCA,$*) \
			/mnt/cwd/$(call FIT_UMAP,$*) \
			/mnt/cwd/quinn_$*.fit.latent_flow.total_euclidean_distance.tsv \
			/mnt/cwd/quinn_$*.fit.latent_flow

clean-gexlineage:
	rm -f *.fit.*

clean-waddington-plots:
	rm -f *.waddington.pdf
