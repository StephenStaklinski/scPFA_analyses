.PHONY:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

TREE_FILE := ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex

# Specify fixed parameters across simulations
N_GENES := 100
K := 3
SIGMA_LATENT := 2.0,1.0,0.5

# Specify the range of sigma-obs values to test
SIGMA_OBS_VALUES := 0.01 0.5 1.0 5.0 10.0 20.0

# Specify how many sims to run for each condition
N_SIMS := 10

FIRST_SIGMA := $(firstword $(SIGMA_OBS_VALUES))
FIRST_EVAL := sim_$(FIRST_SIGMA)_eval.eval.summary.tsv
EVALS := $(addprefix sim_,$(addsuffix _eval.eval.summary.tsv, $(SIGMA_OBS_VALUES)))

all: diff.txt

# Simulate the data under the specified conditions
sim_%.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_$* \
		--k ${K} \
		--n-genes ${N_GENES} \
		--sigma-obs $* \
		--sigma-latent ${SIGMA_LATENT} > sim_$*.term

# Fit the model to the simulated data
sim_%.fit.model.summary.tsv: ${TREE_FILE} sim_%.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_$*.expr.tsv \
		--outprefix sim_$*_fit \
		--no-filter > sim_$*_fit.term

# Evaluate the model recovery of the simulated data
sim_%_eval.eval.summary.tsv: sim_%.truth.summary.tsv sim_%.fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_$* \
		--fit-prefix sim_$*_fit \
		--outprefix sim_$*_eval > sim_$*_eval.term

# Compare recovery metrics
diff.txt: $(EVALS)
	printf "condition" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s", $$1 } END { printf "\n" }' $(FIRST_EVAL) >> $@; \
	for condition in $(SIGMA_OBS_VALUES); do \
		printf "$$condition" >> $@; \
		awk -F '\t' 'NR > 1 { printf "\t%s", $$2 } END { printf "\n" }' sim_$${condition}_eval.eval.summary.tsv >> $@; \
	done

clean:
	rm -f sim_* diff.txt
