.PHONY:

GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage
INPUT_DATA_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets

RELATIVE_DATA_DIR := ${INPUT_DATA_DIR}/jiawei_phylo_signal_gene_filter_sim_data_examples_4_1_26

JIAWEI_POS_RESULTS := ${RELATIVE_DATA_DIR}/plasticity_A_nb.tsv
JIAWEI_NEG_RESULTS := ${RELATIVE_DATA_DIR}/plasticity_A_shuffled_nb.tsv

all: diff.txt

# Test positive examples
positive.correlation.lrt.tsv positive.time: ${RELATIVE_DATA_DIR}/tree.nex ${RELATIVE_DATA_DIR}/readcounts_A.tsv
	/usr/bin/time -o positive.lrt.time /usr/bin/time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${RELATIVE_DATA_DIR}/tree.nex \
		--expr ${RELATIVE_DATA_DIR}/readcounts_A.tsv \
		--outprefix positive \
		--filter-only > positive.term

# Test negative examples
negative.correlation.lrt.tsv negative.time: ${RELATIVE_DATA_DIR}/tree.nex ${RELATIVE_DATA_DIR}/readcounts_A_shuffled.tsv
	/usr/bin/time -o negative.lrt.time /usr/bin/time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees ${RELATIVE_DATA_DIR}/tree.nex \
		--expr ${RELATIVE_DATA_DIR}/readcounts_A_shuffled.tsv \
		--outprefix negative \
		--filter-only > negative.term

diff.txt: positive.correlation.lrt.tsv negative.correlation.lrt.tsv
	total_sims=$$(tail -n +2 positive.correlation.lrt.tsv | wc -l); \
	jiawei_pos_tp=$$(tail -n +2 ${JIAWEI_POS_RESULTS} | cut -f15 | grep -c '^True$$'); \
	jiawei_neg_tn=$$(tail -n +2 ${JIAWEI_NEG_RESULTS} | cut -f15 | grep -c '^False$$'); \
	gexLineage_pos_tp=$$(tail -n +2 positive.correlation.lrt.tsv | cut -f8 | grep -c '^yes$$'); \
	gexLineage_neg_tn=$$(tail -n +2 negative.correlation.lrt.tsv | cut -f8 | grep -c '^no$$'); \
	echo "total_sims,jiawei_pos_tp,jiawei_neg_tn,gexLineage_pos_tp,gexLineage_neg_tn" > $@; \
	echo "$$total_sims,$$jiawei_pos_tp,$$jiawei_neg_tn,$$gexLineage_pos_tp,$$gexLineage_neg_tn" >> $@

    
clean:
	rm -f positive.* negative.*
