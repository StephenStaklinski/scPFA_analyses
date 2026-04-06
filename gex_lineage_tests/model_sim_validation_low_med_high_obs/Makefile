.PHONY:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

TREE_FILE := ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex

N_GENES := 100

all: diff.txt

# Low-noise latent simulation and recovery
sim_test_low.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_test_low \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 0.01 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 7

sim_test_low_fit.model.summary.tsv: ${TREE_FILE} sim_test_low.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_test_low.expr.tsv \
		--outprefix sim_test_low_fit \
		--no-filter \
		--seed 7 > sim_test_low_fit.term

sim_test_low_eval.eval.summary.tsv: sim_test_low.truth.summary.tsv sim_test_low_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_test_low \
		--fit-prefix sim_test_low_fit \
		--outprefix sim_test_low_eval

# Moderate-noise latent simulation and recovery
sim_test_med.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_test_med \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 0.5 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 9

sim_test_med_fit.model.summary.tsv: ${TREE_FILE} sim_test_med.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_test_med.expr.tsv \
		--outprefix sim_test_med_fit \
		--no-filter \
		--seed 9 > sim_test_med_fit.term

sim_test_med_eval.eval.summary.tsv: sim_test_med.truth.summary.tsv sim_test_med_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_test_med \
		--fit-prefix sim_test_med_fit \
		--outprefix sim_test_med_eval

# High-noise latent simulation and recovery
sim_test_high.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_test_high \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 1.0 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 11

sim_test_high_fit.model.summary.tsv: ${TREE_FILE} sim_test_high.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_test_high.expr.tsv \
		--outprefix sim_test_high_fit \
		--no-filter \
		--seed 11 > sim_test_high_fit.term

sim_test_high_eval.eval.summary.tsv: sim_test_high.truth.summary.tsv sim_test_high_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_test_high \
		--fit-prefix sim_test_high_fit \
		--outprefix sim_test_high_eval

# Compare recovery metrics
diff.txt: sim_test_low_eval.eval.summary.tsv sim_test_med_eval.eval.summary.tsv sim_test_high_eval.eval.summary.tsv
	printf "condition" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s", $$1 } END { printf "\n" }' sim_test_low_eval.eval.summary.tsv >> $@; \
	for condition in low med high; do \
		printf "$$condition" >> $@; \
		awk -F '\t' 'NR > 1 { printf "\t%s", $$2 } END { printf "\n" }' sim_test_$${condition}_eval.eval.summary.tsv >> $@; \
	done

clean:
	rm -f sim_test_low* sim_test_med* sim_test_high* diff.txt
