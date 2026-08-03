export SHELL=/usr/bin/bash

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := ../data

PREV_OUTPREFIX := ./prev_run/quinn.90

# Selected the smallest CP with perfectly matching tree and expr matrix cell names
OUTPREFIX := quinn.90
FILTER_PREFIX := ${OUTPREFIX}.filter
FILTERED_GEX := ${FILTER_PREFIX}.filtered.tsv
PCA_INIT := ${OUTPREFIX}.pca.eigenvectors.tsv

all: run compare
run: ${OUTPREFIX}.log
compare: diff-objective.txt diff-runtime.txt

# Subset to less genes for a faster run
subset_${OUTPREFIX}.gex.tsv: ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.gex.tsv
	cut -f1-21 $< > $@

# Filter genes once into the modeling-ready matrix consumed by gexFactor.
${FILTERED_GEX}: ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex subset_${OUTPREFIX}.gex.tsv
	/usr/bin/time -o ${FILTER_PREFIX}.time \
		${GEX_LINEAGE_DIR}/bin/gexFilter \
		--trees $< \
		--expr subset_${OUTPREFIX}.gex.tsv \
		--outprefix ${FILTER_PREFIX} \
		--filter-test lrt-lambda \
		--seed 1 > ${FILTER_PREFIX}.term

# Compute the maxPhyloPCA initialization separately from model fitting.
${PCA_INIT}: ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex ${FILTERED_GEX}
	${GEX_LINEAGE_DIR}/bin/gexPca \
		--tree $< \
		--expr ${FILTERED_GEX} \
		--outprefix ${OUTPREFIX} \
		--dim 5 > ${OUTPREFIX}.pca.term

# Fit the factor model to the filtered matrix.
${OUTPREFIX}.term ${OUTPREFIX}.time ${OUTPREFIX}.log ${OUTPREFIX}.summary.tsv \
${OUTPREFIX}.X.tsv ${OUTPREFIX}.F.tsv ${OUTPREFIX}.L.tsv: \
${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex ${FILTERED_GEX} ${PCA_INIT}
	/usr/bin/time -o ${OUTPREFIX}.time \
		${GEX_LINEAGE_DIR}/bin/gexFactor \
		--trees ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex \
		--expr ${FILTERED_GEX} \
		--outprefix ${OUTPREFIX} \
		--dim 5 \
		--pca ${PCA_INIT} \
		--seed 1 > ${OUTPREFIX}.term

# Compare the final objective value
diff-objective.txt: ${OUTPREFIX}.log ${PREV_OUTPREFIX}.log 
	prev_obj=$$(awk -F '\t' '$$1 == "# best_state" { obj = $$3 } END { print obj }' ${PREV_OUTPREFIX}.log); \
	obj=$$(awk -F '\t' '$$1 == "# best_state" { obj = $$3 } END { print obj }' ${OUTPREFIX}.log); \
	printf "current_objective\tprevious_objective\trelative_diff\n" > $@; \
	awk -v obj="$$obj" -v prev_obj="$$prev_obj" 'BEGIN { rel = (obj - prev_obj) / prev_obj; if (rel == 0) rel = 0; printf "%.17g\t%.17g\t%.6g\n", obj, prev_obj, rel }' >> $@

# Compare runtime
diff-runtime.txt: ${OUTPREFIX}.time ${PREV_OUTPREFIX}.time
	prev_time=$$(head -n 1 ${PREV_OUTPREFIX}.time | cut -d' ' -f1 | sed 's/user//'); \
	time=$$(head -n 1 ${OUTPREFIX}.time | cut -d' ' -f1 | sed 's/user//'); \
	echo -e "current_time\tprevious_time\trelative_diff" > $@; \
	echo -e "$$time\t$$prev_time\t$$(echo "scale=6; ($$time - $$prev_time) / $$prev_time" | bc)" >> $@

# Reset the stored previous run to the current run for future comparisons, and record date/time of the stored run
store-curr-as-prev:
	mkdir -p prev_run
	mv ${OUTPREFIX}.* prev_run/
	date > prev_run/run_stored_date.txt

# Remove the output files from the current run
clean:
	rm -f ${OUTPREFIX}.* diff*.txt
