
MAIN_DIR=/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/gene_phylo_filter_validation_sims

# Get all performance results
OUTFILE=$MAIN_DIR/eval.all.performance.txt
rm -f $OUTFILE

for file in $MAIN_DIR/*/eval.all.performance.txt; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done

# Get all time results
OUTFILE=$MAIN_DIR/eval.all.time.txt
rm -f $OUTFILE

for file in $MAIN_DIR/*/eval.all.time.txt; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done

# Get all stats results
OUTFILE=$MAIN_DIR/eval.all.tree_stats.tsv
rm -f $OUTFILE

for file in $MAIN_DIR/*/eval.all.tree_stats.tsv; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done
