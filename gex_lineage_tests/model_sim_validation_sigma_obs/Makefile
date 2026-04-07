.PHONY:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

TREE_FILE := ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex

# Specify fixed parameters across simulations
N_GENES := 100
K := 3	# Number of latent factors (K) for the simulation
SIGMA_LATENT := 2.0,1.0,0.5	# Comma-sep list of doubles of length K for the latent Brownian process noise

# Specify the range of sigma-obs values to test
SIGMA_OBS_VALUES := 0.01 # 0.5 1.0 5.0 10.0 20.0 50.0 100.0

# Specify how many sims to run for each condition
N_SIMS := 1

SIM_IDS := $(shell seq 1 $(N_SIMS))
RUN_STEMS := $(foreach condition,$(SIGMA_OBS_VALUES),$(foreach sim,$(SIM_IDS),$(condition)_$(sim)))
FIRST_STEM := $(firstword $(RUN_STEMS))
FIRST_EVAL := sim_$(FIRST_STEM).eval.summary.tsv
EVALS := $(addprefix sim_,$(addsuffix .eval.summary.tsv, $(RUN_STEMS)))

all: diff.txt

# Simulate the data under the specified conditions
sim_%.sim.summary.tsv: ${TREE_FILE}
	stem="$*"; condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_$${condition}_$${sim_num}.sim \
		--k ${K} \
		--n-genes ${N_GENES} \
		--sigma2-obs $${condition} \
		--sigma2-latent ${SIGMA_LATENT} \
		--seed $$(shuf -i 1-1000000000 -n 1) > sim_$${condition}_$${sim_num}.sim.term

# Fit the model to the simulated data
sim_%.fit.model.summary.tsv: ${TREE_FILE} sim_%.sim.summary.tsv
	stem="$*"; condition="$${stem%_*}"; sim_num="$${stem##*_}"; \
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_$${condition}_$${sim_num}.sim.expr.tsv \
		--outprefix sim_$${condition}_$${sim_num}.fit \
		--no-filter \
		--seed $$(shuf -i 1-1000000000 -n 1) > sim_$${condition}_$${sim_num}.fit.term

# Evaluate the model recovery of the simulated data
sim_%.eval.summary.tsv: sim_%.sim.summary.tsv sim_%.fit.model.summary.tsv
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
