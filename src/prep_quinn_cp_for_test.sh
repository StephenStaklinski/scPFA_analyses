
# Need to activate environment "vine" on evolgen to have the dependencies for the scripts used here
CP="90"
VINENWK="/local/storage/no-backup/vine-benchmarks/crispr_real_data/quinn-all/quinn.${CP}.var.nwk"
COLLAPSINGMAP="/local/storage/no-backup/vine-benchmarks/crispr_real_data/quinn-all/quinn.${CP}.indels.collapsing_map.tsv"

check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File not found: $1"
        exit 1
    fi
}

check_file "$VINENWK"
check_file "$COLLAPSINGMAP"

OUTPUTDIR="/home/staklins/projects/gex_lineage_project/gex_lineage_tests/input/"

SCRIPTSDIR="/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/src/"


python ${SCRIPTSDIR}/nwk2nex.py $VINENWK ${OUTPUTDIR}/quinn.${CP}.var.nex
cut -f1-2 ${COLLAPSINGMAP} | cut -d',' -f1 > ${OUTPUTDIR}/quinn.${CP}.indels.collapsing_map.dedup.tsv    # Chooses the first barcode of a collapsed set (for now this is okay, but may need to fix this in the future)
python ${SCRIPTSDIR}/rename_nexus_tips.py ${OUTPUTDIR}/quinn.${CP}.var.nex ${OUTPUTDIR}/quinn.${CP}.indels.collapsing_map.dedup.tsv ${OUTPUTDIR}/quinn.${CP}.var.decollapsed.nex
mv ${OUTPUTDIR}/quinn.${CP}.var.decollapsed.nex ${OUTPUTDIR}/quinn.${CP}.var.nex
cut -f2 ${OUTPUTDIR}/quinn.${CP}.indels.collapsing_map.dedup.tsv > ${OUTPUTDIR}/quinn.${CP}.barcodes.txt
rm ${OUTPUTDIR}/quinn.${CP}.indels.collapsing_map.dedup.tsv
python ${SCRIPTSDIR}/subset_cells_10x_mtx_to_csv.py \
    /home/staklins/projects/gex_lineage_project/gex_lineage_tests/input/GSM4905335_barcodes.5k.tsv \
    /home/staklins/projects/gex_lineage_project/gex_lineage_tests/input/GSM4905335_genes.5k.tsv \
    /home/staklins/projects/gex_lineage_project/gex_lineage_tests/input/GSM4905335_matrix.5k.mtx \
    ${OUTPUTDIR}/quinn.${CP}.barcodes.txt \
    ${OUTPUTDIR}/quinn.${CP}.gex.tsv
rm ${OUTPUTDIR}/quinn.${CP}.barcodes.txt


# Maybe need to add a rule to remove missing cells from the tree as well (could just do it in the C code..., but first maybe need to figure out why cells are missing at all)

