.PHONY: 

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

PREV_OUTPREFIX := ./prev_run/quinn.90

# Selected the smallest CP with perfectly matching tree and expr matrix cell names
OUTPREFIX := quinn.90

# Quinn mouse 5k experiment time
TREE_RESCALE_TOTAL_TIME := 54

all: run compare
run: ${OUTPREFIX}.model.log
compare: diff-objective.txt diff-runtime.txt

# Main run on the example dataset
${OUTPREFIX}.term ${OUTPREFIX}.time ${OUTPREFIX}.model.log ${OUTPREFIX}.model.summary.tsv \
${OUTPREFIX}.model.Z.tsv ${OUTPREFIX}.model.L.tsv: \
${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.gex.tsv
	/usr/bin/time -o ${OUTPREFIX}.time \
		${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.var.nex \
		--expr ${INPUT_DATA_DIR}/quinn_CP90/${OUTPREFIX}.gex.tsv \
		--outprefix ${OUTPREFIX} \
		--seed 1 \
		--tree-total-time ${TREE_RESCALE_TOTAL_TIME} > ${OUTPREFIX}.term

# Compare the final objective value
diff-objective.txt: ${OUTPREFIX}.model.log ${PREV_OUTPREFIX}.model.log 
	prev_obj=$$(tail -n 2 ${PREV_OUTPREFIX}.model.log | head -n 1 | cut -f2); \
	obj=$$(tail -n 2 ${OUTPREFIX}.model.log | head -n 1 | cut -f2); \
	echo -e "current_objective\tprevious_objective\trelative_diff" > $@; \
	echo -e "$$obj\t$$prev_obj\t$$(echo "scale=6; ($$obj - $$prev_obj) / $$prev_obj" | bc)" >> $@

# Compare runtime
diff-runtime.txt: ${OUTPREFIX}.time ${PREV_OUTPREFIX}.time
	prev_time=$$(head -n 1 ${PREV_OUTPREFIX}.time | cut -d' ' -f1 | sed 's/user//'); \
	time=$$(head -n 1 ${OUTPREFIX}.time | cut -d' ' -f1 | sed 's/user//'); \
	echo -e "current_time\tprevious_time\trelative_diff" > $@; \
	echo -e "$$time\t$$prev_time\t$$(echo "scale=6; ($$time - $$prev_time) / $$prev_time" | bc)" >> $@

# Reset the stored previous run to the current run for future comparisons
store-curr-as-prev:
	mkdir -p prev_run
	mv ${OUTPREFIX}.* prev_run/

# Remove the output files from the current run
clean:
	rm -f ${OUTPREFIX}.*