.PHONY:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

TREE_FILE := ${INPUT_DATA_DIR}/quinn_CP90/quinn.90.var.nex

N_GENES := 100

all: diff.txt

# Low-noise latent simulation and recovery
sim_low.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_low \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 0.01 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 7

sim_low_fit.model.summary.tsv: ${TREE_FILE} sim_low.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_low.expr.tsv \
		--outprefix sim_low_fit \
		--no-filter \
		--seed 7 > sim_low_fit.term

sim_low_eval.eval.summary.tsv: sim_low.truth.summary.tsv sim_low_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_low \
		--fit-prefix sim_low_fit \
		--outprefix sim_low_eval

# Moderate-noise latent simulation and recovery
sim_med.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_med \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 1.0 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 9

sim_med_fit.model.summary.tsv: ${TREE_FILE} sim_med.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_med.expr.tsv \
		--outprefix sim_med_fit \
		--no-filter \
		--seed 9 > sim_med_fit.term

sim_med_eval.eval.summary.tsv: sim_med.truth.summary.tsv sim_med_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_med \
		--fit-prefix sim_med_fit \
		--outprefix sim_med_eval

# High-noise latent simulation and recovery
sim_high.truth.summary.tsv: ${TREE_FILE}
	${GEX_LINEAGE_DIR}/bin/gexSimulateLatent \
		--trees ${TREE_FILE} \
		--outprefix sim_high \
		--k 3 \
		--n-genes ${N_GENES} \
		--sigma-obs 5.0 \
		--sigma-latent 2.0,1.0,0.5 \
		--seed 11

sim_high_fit.model.summary.tsv: ${TREE_FILE} sim_high.truth.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${TREE_FILE} \
		--expr sim_high.expr.tsv \
		--outprefix sim_high_fit \
		--no-filter \
		--seed 11 > sim_high_fit.term

sim_high_eval.eval.summary.tsv: sim_high.truth.summary.tsv sim_high_fit.model.summary.tsv
	${GEX_LINEAGE_DIR}/bin/gexEvaluateSimulation \
		--truth-prefix sim_high \
		--fit-prefix sim_high_fit \
		--outprefix sim_high_eval

# Compare recovery metrics
diff.txt: sim_low_eval.eval.summary.tsv sim_med_eval.eval.summary.tsv sim_high_eval.eval.summary.tsv
	printf "condition" > $@; \
	awk -F '\t' 'NR > 1 { printf "\t%s", $$1 } END { printf "\n" }' sim_low_eval.eval.summary.tsv >> $@; \
	for condition in low med high; do \
		printf "$$condition" >> $@; \
		awk -F '\t' 'NR > 1 { printf "\t%s", $$2 } END { printf "\n" }' sim_$${condition}_eval.eval.summary.tsv >> $@; \
	done

clean:
	rm -f sim_* diff.txt
