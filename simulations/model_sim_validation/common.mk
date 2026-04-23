export SHELL=/usr/bin/bash

.SECONDARY:
.PRECIOUS:

MAIN_DIR := /home/staklins/projects/gex_lineage_project
BENCHMARKS_DIR := $(MAIN_DIR)/gex_lineage_benchmarks
REL_PATH := simulations/model_sim_validation/$(NTAXA)taxa$(NGENES)genes
CURR_DIR := $(BENCHMARKS_DIR)/$(REL_PATH)

CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := ${CONTAINERS}/pyplotting/pyplotting.sif
CASSIOPEIA_SIF := ${CONTAINERS}/cassiopeia/cassiopeia.sif

SRC := ${BENCHMARKS_DIR}/src
GEX_LINEAGE_DIR := $(MAIN_DIR)/gex_lineage


SIM_IDS := $(shell seq 1 $(NSAMP))
TREES := $(shell seq -f sim.%.0f.tree.nex 1 $(NSAMP))
SIMS := $(shell seq -f sim.%.0f.sim.summary.tsv 1 $(NSAMP))
FITS := $(shell seq -f sim.%.0f.fit.summary.tsv 1 $(NSAMP))
TIMES := $(shell seq -f sim.%.0f.fit.time 1 $(NSAMP))
PANELFIGS := $(shell seq -f sim.%.0f.panel.png 1 $(NSAMP))
ALIGNEDFFIGS := $(shell seq -f sim.%.0f.alignedF.png 1 $(NSAMP))
FITS_LOG_PDFS := $(shell seq -f sim.%.0f.fit.log.pdf 1 $(NSAMP))
FIRST_EVAL := sim.1.eval.summary.tsv
EVALS := $(shell seq -f sim.%.0f.eval.summary.tsv 1 $(NSAMP))

all: simulate fit eval
simulate: $(TREES) $(SIMS)
fit: $(FITS) $(FITS_LOG_PDFS) eval.all.fit.summaries.pdf $(PANELFIGS) $(ALIGNEDFFIGS)
eval: $(EVALS) diff.txt eval.all.time.txt

# Simulate the cell lineage tree
sim.%.tree.nwk:
	singularity exec --bind $(BENCHMARKS_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/simulateCellTree.py \
		--out_tree /mnt/$(REL_PATH)/sim.$*.tree.nwk \
		--num_tips $(NTAXA) \
		--birth_rate 0.075 \
		--death_rate 0.005 \
		--desired_time $(TOTAL_TIME)

# Convert the simulated tree from Newick to Nexus format
sim.%.tree.nex: sim.%.tree.nwk
	singularity exec --bind $(BENCHMARKS_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py \
		/mnt/$(REL_PATH)/sim.$*.tree.nwk \
		/mnt/$(REL_PATH)/sim.$*.tree.nex

sim.%.sim.summary.tsv sim.%.sim.F.tsv sim.%.sim.L.tsv sim.%.sim.X.tsv: sim.%.tree.nex
	${GEX_LINEAGE_DIR}/bin/gexSim \
		--seed $$(shuf -i 1-1000000000 -n 1) \
		--trees sim.$*.tree.nex \
		--outprefix sim.$*.sim \
		--tree-total-time $(TOTAL_TIME) \
		--n-genes ${NGENES} \
		--L-l2-norm ${DESIRED_L_ROW_NORMS} \
		--use-n-trees -1 \
		--dim ${K} \
		--sigma2-obs $(SIGMA2_OBS) \
		--include-factorization > sim.$*.sim.term

# Fit the model to the simulated data
sim.%.fit.time sim.%.fit.summary.tsv sim.%.fit.log sim.%.fit.F.tsv sim.%.fit.L.tsv sim.%.fit.X.tsv: sim.%.tree.nex sim.%.sim.summary.tsv
	/usr/bin/time -o sim.$*.fit.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--seed $$(shuf -i 1-1000000000 -n 1) \
		--trees sim.$*.tree.nex \
		--expr sim.$*.sim.X.tsv \
		--outprefix sim.$*.fit \
		--dim ${K} \
		--no-filter \
		--no-preprocess \
		--seed $$(shuf -i 1-1000000000 -n 1) > sim.$*.fit.term

# Plot the fit optimization log results
sim.%.fit.log.pdf: sim.%.fit.log
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_optimization_log.py /mnt/cwd/sim.$*.fit.log -1 0 /mnt/cwd/sim.$*.fit.log

sim.%.panel.png: sim.%.sim.F.tsv sim.%.sim.L.tsv sim.%.sim.X.tsv sim.%.fit.F.tsv sim.%.fit.L.tsv sim.%.fit.X.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_sim_fit_panel.py \
			/mnt/cwd/sim.$*.sim.F.tsv \
			/mnt/cwd/sim.$*.sim.L.tsv \
			/mnt/cwd/sim.$*.sim.X.tsv \
			/mnt/cwd/sim.$*.fit.F.tsv \
			/mnt/cwd/sim.$*.fit.L.tsv \
			/mnt/cwd/sim.$*.fit.X.tsv \
			/mnt/cwd/sim.$*.panel.png

sim.%.alignedF.png: sim.%.sim.F.tsv sim.%.fit.F.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/procrustes_align_plot_F.py \
			/mnt/cwd/sim.$*.sim.F.tsv \
			/mnt/cwd/sim.$*.fit.F.tsv \
			/mnt/cwd/sim.$*.alignedF.png

eval.all.fit.summaries.tsv: $(FITS) $(SIMS)
	printf "sim_num" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s_fit\t%s_simulated", $$1, $$1 } END { printf "\n" }' sim.1.fit.summary.tsv >> $@; \
	for sim_num in $(SIM_IDS); do \
		printf "%s" "$$sim_num" >> $@; \
		paste sim.$${sim_num}.fit.summary.tsv sim.$${sim_num}.sim.summary.tsv | \
		awk -F '\t' 'NR > 1 { printf "\t%s\t%s", $$2, $$4 } END { printf "\n" }' >> $@; \
	done

eval.all.fit.summaries.pdf: eval.all.fit.summaries.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_sim_fit_summary_reps.py /mnt/cwd/eval.all.fit.summaries.tsv /mnt/cwd/eval.all.fit.summaries.pdf

sim.%.eval.summary.tsv: sim.%.sim.summary.tsv sim.%.fit.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvalSim \
		--sim-prefix sim.$*.sim \
		--fit-prefix sim.$*.fit \
		--outprefix sim.$*.eval > sim.$*.eval.term

diff.txt: $(EVALS)
	printf "NTAXA\tNGENES\tsim" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s", $$1 } END { printf "\n" }' $(FIRST_EVAL) >> $@; \
	for sim_num in $(SIM_IDS); do \
		printf "%s\t%s\t%s" "$(NTAXA)" "$(NGENES)" "$$sim_num" >> $@; \
		awk -F '\t' 'NR > 1 { printf "\t%s", $$2 } END { printf "\n" }' sim.$${sim_num}.eval.summary.tsv >> $@; \
	done

eval.all.time.txt: $(TIMES)
	printf "NTAXA\tNGENES\tsim_num\tsec\n" > $@; \
	for sim_num in $(SIM_IDS); do \
		fit_time=$$(head -n 1 "sim.$${sim_num}.fit.time" | cut -d" " -f1 | sed 's/user//g'); \
		printf "%s\t%s\t%s\t%s\n" "$(NTAXA)" "$(NGENES)" "$$sim_num" "$$fit_time" >> $@; \
	done

clean:
	rm -f sim.* diff.txt eval.all*

archive-all:
	archive_dir=archive_all_$(shell date +%Y-%m-%d_%H-%M); \
	mkdir -p $$archive_dir; \
	mv sim* $$archive_dir/; \
	mv eval.all* $$archive_dir/; \
	mv diff.txt $$archive_dir/;

archive-fit:
	archive_dir=archive_fit_$(shell date +%Y-%m-%d_%H-%M); \
	mkdir -p $$archive_dir; \
	mv *.fit* $$archive_dir/; \
	mv *.alignedF.png $$archive_dir/; \
	mv *.panel.png $$archive_dir/; \
	mv *.eval* $$archive_dir/; \
	mv eval.* $$archive_dir/; \
	mv diff.txt $$archive_dir/;
