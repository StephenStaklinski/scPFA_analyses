
MAIN_DIR=/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/simulations/model_sim_validation


# Get diff files from gexEvalSim
OUTFILE=$MAIN_DIR/diff.txt
rm -f $OUTFILE

diff_files=$(find $MAIN_DIR -mindepth 2 -maxdepth 2 ! -path "$MAIN_DIR/archive*" -type f -name "diff.txt")
for file in $diff_files; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done

python $MAIN_DIR/plot_diff.py $OUTFILE $MAIN_DIR/diff.pdf

# Get sim and fit summaries
OUTFILE=$MAIN_DIR/eval.all.fit.summaries.tsv
rm -f $OUTFILE

summary_files=$(find $MAIN_DIR -mindepth 2 -maxdepth 2 ! -path "$MAIN_DIR/archive*" -type f -name "eval.all.fit.summaries.tsv")
for file in $summary_files; do
    if [[ ! -f $OUTFILE ]]; then
        # Write the header from the first file
        head -n 1 $file > $OUTFILE
    fi

    # Append the data from the current file, skipping the header
    tail -n +2 $file >> $OUTFILE
done

python $MAIN_DIR/plot_summaries.py $OUTFILE $MAIN_DIR/summaries.pdf


# Get runtimes
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

python $MAIN_DIR/plot_all_time.py $OUTFILE $MAIN_DIR/times.pdf
