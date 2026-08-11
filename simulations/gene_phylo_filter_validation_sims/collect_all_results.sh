
MAIN_DIR=/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/simulations/gene_phylo_filter_validation_sims

# Get all performance results
OUTFILE=$MAIN_DIR/eval.all.performance.txt
rm -f $OUTFILE

performance_files=$(find $MAIN_DIR -mindepth 2 -maxdepth 2 ! -path "$MAIN_DIR/archive*" -type f -name "eval.all.performance.txt")
for file in $performance_files; do
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

time_files=$(find $MAIN_DIR -mindepth 2 -maxdepth 2 ! -path "$MAIN_DIR/archive*" -type f -name "eval.all.time.txt")
for file in $time_files; do
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

tree_stats_files=$(find $MAIN_DIR -mindepth 2 -maxdepth 2 ! -path "$MAIN_DIR/archive*" -type f -name "eval.all.tree_stats.tsv")
for file in $tree_stats_files; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done
