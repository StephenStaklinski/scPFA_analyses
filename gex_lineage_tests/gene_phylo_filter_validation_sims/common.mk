export SHELL=/usr/bin/bash

.PRECIOUS: tree.%.true.nwk tree.%.correlation.path.tsv


GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage

# Absolute path
MAIN_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks

REL_PATH := gex_lineage_tests/gene_phylo_filter_validation_sims

# Relative paths
CONTAINERS := $(MAIN_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
PATH_SIF := $(CONTAINERS)/path/path.sif

TREES := $(shell seq -f tree.%.0f.true.nex 1 $(NSAMP))
EXPRS := $(shell seq -f tree.%.0f.true.expr.tsv 1 $(NSAMP))
LAMBDATERMS := $(shell seq -f tree.%.0f.correlation.lrt.lambda.term 1 $(NSAMP))
LAMBDATIMES := $(shell seq -f tree.%.0f.correlation.lrt.lambda.time 1 $(NSAMP))
FULLTERMS := $(shell seq -f tree.%.0f.correlation.lrt.full.term 1 $(NSAMP))
FULLTIMES := $(shell seq -f tree.%.0f.correlation.lrt.full.time 1 $(NSAMP))
MORANTERMS := $(shell seq -f tree.%.0f.correlation.moran.term 1 $(NSAMP))
MORANTIMES := $(shell seq -f tree.%.0f.correlation.moran.time 1 $(NSAMP))
PATHTERMS := $(shell seq -f tree.%.0f.correlation.path.term 1 $(NSAMP))
PATHTIMES := $(shell seq -f tree.%.0f.correlation.path.time 1 $(NSAMP))

NSIMS := 100

all: eval.all.performance.txt eval.all.time.txt

simulate: $(TREES) $(EXPRS)
lambda: $(LAMBDATERMS)
full: $(FULLTERMS)
moran: $(MORANTERMS)
path: $(PATHTERMS)

# Simulate the cell lineage trees
tree.%.true.nwk:
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/simulateCellTree.py \
		--out_tree /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nwk \
		--num_tips $(NTAXA) \
		--birth_rate 0.075 \
		--death_rate 0.005 \
		--desired_time 54

# Convert the simulated trees from Newick to Nexus format
tree.%.true.nex: tree.%.true.nwk
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nwk \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nex

# Generate a dummy gene expression matrix (required by gexLineage, though it is not used)
tree.%.true.expr.tsv: tree.%.true.nex
	awk 'BEGIN { inblock=0; print "cell\tdummy" } /TaxLabels/ { inblock=1; next } inblock { line=$$0; gsub(/^[ \t]+|[ \t]+$$/, "", line); has_end = (line ~ /;$$/); sub(/;$$/, "", line); if (line != "") print line "\t0"; if (has_end) inblock=0; }' $< > $@

# Run the Pagel's lambda LRT phylo signal filter
tree.%.correlation.lrt.lambda.tsv tree.%.correlation.lrt.lambda.term tree.%.correlation.lrt.lambda.time tree.%.phylo_filter_sims.expr.tsv: tree.%.true.nex tree.%.true.expr.tsv
	/usr/bin/time -o tree.$*.correlation.lrt.lambda.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.expr.tsv \
		--outprefix tree.$* \
		--filter-test lrt \
		--lrt-alt lambda \
		--sim-filter-only \
		--write-filter-sims \
		--n-sims $(NSIMS) \
		--seed $$(shuf -i 1-1000000000 -n 1) > tree.$*.correlation.lrt.lambda.term

# Run the full LRT phylo signal filter
tree.%.correlation.lrt.full.tsv tree.%.correlation.lrt.full.term tree.%.correlation.lrt.full.time: tree.%.true.nex tree.%.true.expr.tsv
	/usr/bin/time -o tree.$*.correlation.lrt.full.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.expr.tsv \
		--outprefix tree.$* \
		--filter-test lrt \
		--lrt-alt full \
		--sim-filter-only \
		--n-sims $(NSIMS) \
		--seed $$(shuf -i 1-1000000000 -n 1) > tree.$*.correlation.lrt.full.term

# Run my implementation of the PATH-based Moran's I autocorrelation phylo signal filter
tree.%.correlation.moran.tsv tree.%.correlation.moran.term tree.%.correlation.moran.time: tree.%.true.nex tree.%.true.expr.tsv
	/usr/bin/time -o tree.$*.correlation.moran.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.expr.tsv \
		--outprefix tree.$* \
		--filter-test moran \
		--sim-filter-only \
		--n-sims $(NSIMS) \
		--seed $$(shuf -i 1-1000000000 -n 1) > tree.$*.correlation.moran.term

# Run PATH's implementation of the autocorrelation phylo signal filter
tree.%.correlation.path.tsv tree.%.correlation.path.time: tree.%.phylo_filter_sims.expr.tsv
	singularity exec --bind $(MAIN_DIR):/mnt $(PATH_SIF) /usr/bin/time -o /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.correlation.path.time \
	Rscript /mnt/src/path_gene_phylo_correlation.R \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nex \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.phylo_filter_sims.expr.tsv \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.correlation.path.tsv

# Summarize PATH results in the same format as the other methods
tree.%.correlation.path.term: tree.%.correlation.path.tsv
	awk -F '\t' 'BEGIN { tp=0; fn=0; tn=0; fp=0 } \
		NR > 1 { \
			if ($$1 ~ /^gene_/) { \
				if ($$6 == "TRUE") tp++; \
				else fn++; \
			} \
			else if ($$1 ~ /^neg_/) { \
				if ($$6 == "TRUE") fp++; \
				else tn++; \
			} \
		} \
		END { \
			printf "  positives simulated: %d, detected: %d, missed: %d\n", tp + fn, tp, fn; \
			printf "  negatives simulated: %d, rejected: %d, false positives: %d\n", tn + fp, tn, fp; \
		}' $< > $@

# Gather all results to one performance summary file
eval.all.performance.txt: $(LAMBDATERMS) $(FULLTERMS) $(MORANTERMS) $(PATHTERMS)
	{ \
		printf "ntaxa\tmethod\tsim_num\tTP\tFN\tTN\tFP\n"; \
		for f in $(LAMBDATERMS) $(FULLTERMS) $(MORANTERMS) $(PATHTERMS); do \
			ntaxa=$(NTAXA); \
			method=$$(basename "$$f" | cut -d"." -f4- | sed 's/\.term//' | tr '.' '_'); \
			sim_num=$$(basename "$$f" | cut -d"." -f2); \
			awk -v ntaxa="$$ntaxa" -v method="$$method" -v sim_num="$$sim_num" '\
				/positives simulated:/ { tp=$$5; gsub(/,/, "", tp); fn=$$7 } \
				/negatives simulated:/ { tn=$$5; gsub(/,/, "", tn); fp=$$8 } \
				END { print ntaxa "\t" method "\t" sim_num "\t" tp "\t" fn "\t" tn "\t" fp }\
			' "$$f"; \
		done; \
	} > $@

# Gather all runtimes to one summary file
eval.all.time.txt: $(LAMBDATIMES) $(FULLTIMES) $(MORANTIMES) $(PATHTIMES)
	{ \
		printf "ntaxa\tmethod\tsim_num\ttime_sec\n"; \
		for f in $(LAMBDATIMES) $(FULLTIMES) $(MORANTIMES) $(PATHTIMES); do \
			ntaxa=$(NTAXA); \
			method=$$(basename "$$f" | cut -d"." -f4- | sed 's/\.time//' | tr '.' '_'); \
			sim_num=$$(basename "$$f" | cut -d"." -f2); \
			time_sec=$$(head -n 1 "$$f" | cut -d" " -f1 | sed 's/user//g'); \
			echo -e "$$ntaxa\t$$method\t$$sim_num\t$$time_sec"; \
		done; \
	} > $@

clean:
	rm -f tree.* eval.all*

