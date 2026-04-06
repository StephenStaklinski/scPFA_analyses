
# Test correlation metrics vs Jiawei's results
# Run from within gex_lineage directory
/usr/bin/time ./bin/gexLineage \
    --trees ../gex_lineage_tests/input/jiawei_corr_test_example/tree.nex \
    --expr ../gex_lineage_tests/input/jiawei_corr_test_example/readcounts_A_shuffled.tsv \
    --outprefix ../gex_lineage_tests/output/jiawei_corr_test_example/result_shuffled \
    --filter-only
