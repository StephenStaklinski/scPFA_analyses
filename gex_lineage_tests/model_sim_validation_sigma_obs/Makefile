export SHELL=/usr/bin/bash

MAIN_DIR := /home/staklins/projects/gex_lineage_project
BENCHMARKS_DIR := $(MAIN_DIR)/gex_lineage_benchmarks
CURR_DIR := $(BENCHMARKS_DIR)/gex_lineage_tests/model_sim_validation_sigma_obs

CONTAINERS := $(BENCHMARKS_DIR)/containers
PYPLOTTING_SIF := ${CONTAINERS}/pyplotting/pyplotting.sif

SRC := ${BENCHMARKS_DIR}/src

GEX_LINEAGE_DIR := $(MAIN_DIR)/gex_lineage
INPUT_DATA_DIR := $(BENCHMARKS_DIR)/gex_lineage_tests/input_datasets

TREE_FILE := ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex

# Specify fixed parameters across simulations
TOTAL_TIME := 54.0
N_GENES := 1000
DESIRED_TIP_VAR := 20.0,10.0,5.0,2.0,1.0
K := 5	# Number of latent factors for the simulation

# Specify the range of sigma-obs values to test
SIGMA_OBS_VALUES := 0.5 # 0.01 0.5 1.0 5.0 10.0 20.0 50.0 100.0

# Specify how many sims to run for each condition
N_SIMS := 3

SIM_IDS := $(shell seq 1 $(N_SIMS))
RUN_STEMS := $(foreach condition,$(SIGMA_OBS_VALUES),$(foreach sim,$(SIM_IDS),$(condition)_$(sim)))
SIMS := $(addprefix sim_,$(addsuffix .sim.summary.tsv, $(RUN_STEMS)))
FITS := $(addprefix sim_,$(addsuffix .fit.summary.tsv, $(RUN_STEMS)))
FITS_LOG_PDFS := $(addprefix sim_,$(addsuffix .fit.log.pdf, $(RUN_STEMS)))
FIRST_STEM := $(firstword $(RUN_STEMS))
FIRST_EVAL := sim_$(FIRST_STEM).eval.summary.tsv
EVALS := $(addprefix sim_,$(addsuffix .eval.summary.tsv, $(RUN_STEMS)))

all: eval.all.fit.summaries.pdf $(FITS_LOG_PDFS) diff.txt
simulate: $(SIMS)
fit: $(FITS) $(FITS_LOG_PDFS)
eval: $(EVALS)

# Simulate the data under the specified conditions
sim_%.sim.summary.tsv: ${TREE_FILE}
	stem="$*"; condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
	${GEX_LINEAGE_DIR}/bin/gexSim \
		--trees ${TREE_FILE} \
		--outprefix sim_$${condition}_$${sim_num}.sim \
		--tree-total-time $(TOTAL_TIME) \
		--n-genes ${N_GENES} \
		--desired-tip-var ${DESIRED_TIP_VAR} \
		--use-n-trees -1 \
		--dim ${K} \
		--sigma2-obs $${condition} \
		--include-factorization > sim_$${condition}_$${sim_num}.sim.term

# Fit the model to the simulated data
sim_%.fit.summary.tsv sim_%.fit.log: ${TREE_FILE} sim_%.sim.summary.tsv
	stem="$*"; condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
	/usr/bin/time -o sim_$${condition}_$${sim_num}.fit.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_$${condition}_$${sim_num}.sim.X.tsv \
		--outprefix sim_$${condition}_$${sim_num}.fit \
		--dim ${K} \
		--no-filter \
		--no-preprocess \
		--seed $$(shuf -i 1-1000000000 -n 1) > sim_$${condition}_$${sim_num}.fit.term

# Plot the fit optimization log results
sim_%.fit.log.pdf: sim_%.fit.log
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_optimization_log.py /mnt/cwd/sim_$*.fit.log -1 0 /mnt/cwd/sim_$*.fit.log

eval.all.fit.summaries.tsv: $(FITS) $(SIMS)
	printf "condition\tsim_num" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s_fit\t%s_simulated", $$1, $$1 } END { printf "\n" }' sim_$(FIRST_STEM).fit.summary.tsv >> $@; \
	for stem in $(RUN_STEMS); do \
		condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
		printf "%s\t%s" "$$condition" "$$sim_num" >> $@; \
		paste sim_$${stem}.fit.summary.tsv sim_$${stem}.sim.summary.tsv | \
		awk -F '\t' 'NR > 1 { printf "\t%s\t%s", $$2, $$4 } END { printf "\n" }' >> $@; \
	done

eval.all.fit.summaries.pdf: eval.all.fit.summaries.tsv
	singularity exec --bind $(SRC):/mnt/src/ --bind $(CURR_DIR):/mnt/cwd/ $(PYPLOTTING_SIF) \
		python /mnt/src/plot_sim_fit_summary_reps.py /mnt/cwd/eval.all.fit.summaries.tsv /mnt/cwd/eval.all.fit.summaries.pdf

# Evaluate the model recovery of the simulated data
sim_%.eval.summary.tsv: sim_%.sim.summary.tsv sim_%.fit.summary.tsv
	stem="$*"; condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--sim-prefix sim_$${condition}_$${sim_num}.sim \
		--fit-prefix sim_$${condition}_$${sim_num}.fit \
		--outprefix sim_$${condition}_$${sim_num}.eval > sim_$${condition}_$${sim_num}.eval.term

# Compare recovery metrics
diff.txt: $(EVALS)
	printf "condition" > $@; \
	printf "\tsim" >> $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s", $$1 } END { printf "\n" }' $(FIRST_EVAL) >> $@; \
	for stem in $(RUN_STEMS); do \
		condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
		printf "%s\t%s" "$$condition" "$$sim_num" >> $@; \
		awk -F '\t' 'NR > 1 { printf "\t%s", $$2 } END { printf "\n" }' sim_$${stem}.eval.summary.tsv >> $@; \
	done

clean:
	rm -f sim_* diff.txt

archive-all:
	archive_dir=archive_all_$(shell date +%Y-%m-%d_%H-%M); \
	mkdir -p $$archive_dir; \
	mv sim* $$archive_dir/

archive-fit:
	archive_dir=archive_fit_$(shell date +%Y-%m-%d_%H-%M); \
	mkdir -p $$archive_dir; \
	mv *.fit* $$archive_dir/
