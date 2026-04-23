export SHELL=/usr/bin/bash

.PRECIOUS: tree.%.true.nwk tree.%.correlation.path.tsv tree.%.true.pos.expr.tsv tree.%.true.neg.expr.tsv


GEX_LINEAGE_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage

# Absolute path
MAIN_DIR := /home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks

REL_PATH := simulations/gene_phylo_filter_validation_sims

# Relative paths
CONTAINERS := $(MAIN_DIR)/containers
CASSIOPEIA_SIF := $(CONTAINERS)/cassiopeia/cassiopeia.sif
PATH_SIF := $(CONTAINERS)/path/path.sif
BIOPYTHON_SIF := $(CONTAINERS)/biopython/biopython.sif

TREES := $(shell seq -f tree.%.0f.true.nex 1 $(NSAMP))
TREESTATS := $(shell seq -f tree.%.0f.true.stats.tsv 1 $(NSAMP))
EXPRSPOS := $(shell seq -f tree.%.0f.true.pos.expr.tsv 1 $(NSAMP))
EXPRSNEG := $(shell seq -f tree.%.0f.true.neg.expr.tsv 1 $(NSAMP))
LAMBDATSVSPOS := $(shell seq -f tree.%.0f.pos.correlation.lrt.lambda.tsv 1 $(NSAMP))
LAMBDATSVSNEG := $(shell seq -f tree.%.0f.neg.correlation.lrt.lambda.tsv 1 $(NSAMP))
LAMBDATIMESPOS := $(shell seq -f tree.%.0f.pos.correlation.lrt.lambda.time 1 $(NSAMP))
LAMBDATIMESNEG := $(shell seq -f tree.%.0f.neg.correlation.lrt.lambda.time 1 $(NSAMP))
FULLTSVSPOS := $(shell seq -f tree.%.0f.pos.correlation.lrt.full.tsv 1 $(NSAMP))
FULLTIMESPOS := $(shell seq -f tree.%.0f.pos.correlation.lrt.full.time 1 $(NSAMP))
FULLTSVSNEG := $(shell seq -f tree.%.0f.neg.correlation.lrt.full.tsv 1 $(NSAMP))
FULLTIMESNEG := $(shell seq -f tree.%.0f.neg.correlation.lrt.full.time 1 $(NSAMP))
MORANTSVSPOS := $(shell seq -f tree.%.0f.pos.correlation.moran.tsv 1 $(NSAMP))
MORANTIMESPOS := $(shell seq -f tree.%.0f.pos.correlation.moran.time 1 $(NSAMP))
MORANTSVSNEG := $(shell seq -f tree.%.0f.neg.correlation.moran.tsv 1 $(NSAMP))
MORANTIMESNEG := $(shell seq -f tree.%.0f.neg.correlation.moran.time 1 $(NSAMP))
PATHTSVSPOS := $(shell seq -f tree.%.0f.pos.correlation.path.tsv 1 $(NSAMP))
PATHTIMESPOS := $(shell seq -f tree.%.0f.pos.correlation.path.time 1 $(NSAMP))
PATHTSVSNEG := $(shell seq -f tree.%.0f.neg.correlation.path.tsv 1 $(NSAMP))
PATHTIMESNEG := $(shell seq -f tree.%.0f.neg.correlation.path.time 1 $(NSAMP))

NGENES := 1000
DESIRED_TIP_VAR := 0.25

all: eval.all.performance.txt eval.all.time.txt # eval.all.tree_stats.tsv

simulate: $(TREES) $(EXPRSPOS) $(EXPRSNEG)
lambda: $(LAMBDATSVSPOS) $(LAMBDATSVSNEG)
full: $(FULLTSVSPOS) $(FULLTSVSNEG)
moran: $(MORANTSVSPOS) $(MORANTSVSNEG)
path: $(PATHTSVSPOS) $(PATHTSVSNEG)

# Simulate the cell lineage tree
tree.%.true.nwk:
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/simulateCellTree.py \
		--out_tree /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nwk \
		--num_tips $(NTAXA) \
		--birth_rate 0.075 \
		--death_rate 0.005 \
		--desired_time 54

# Get statistics describing the simulated tree
tree.%.true.stats.tsv: tree.%.true.nwk
	singularity exec --bind $(MAIN_DIR):/mnt $(BIOPYTHON_SIF) \
	python /mnt/src/tree_metrics.py /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nwk /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.stats.tsv

# Convert the simulated tree from Newick to Nexus format
tree.%.true.nex: tree.%.true.nwk
	singularity exec --bind $(MAIN_DIR):/mnt $(CASSIOPEIA_SIF) \
	python /mnt/src/nwk2nex.py \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nwk \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nex

# Simulate the positive gene expression matrix
tree.%.true.pos.expr.tsv: tree.%.true.nex
	${GEX_LINEAGE_DIR}/bin/gexSim \
		--trees tree.$*.true.nex \
		--outprefix tree.$*.true.pos \
		--n-genes $(NGENES) \
		--desired-tip-var $(DESIRED_TIP_VAR)

# Simulate the negative gene expression matrix
tree.%.true.neg.expr.tsv: tree.%.true.nex
	${GEX_LINEAGE_DIR}/bin/gexSim \
		--trees tree.$*.true.nex \
		--outprefix tree.$*.true.neg \
		--n-genes $(NGENES) \
		--desired-tip-var $(DESIRED_TIP_VAR) \
		--identity-cov

# Run the Pagel's lambda LRT phylo signal filter on pos sims
tree.%.pos.correlation.lrt.lambda.tsv tree.%.pos.correlation.lrt.lambda.time: tree.%.true.nex tree.%.true.pos.expr.tsv
	/usr/bin/time -o tree.$*.pos.correlation.lrt.lambda.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.pos.expr.tsv \
		--outprefix tree.$*.pos \
		--filter-test lrt \
		--lrt-alt lambda \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run the Pagel's lambda LRT phylo signal filter on neg sims
tree.%.neg.correlation.lrt.lambda.tsv tree.%.neg.correlation.lrt.lambda.time: tree.%.true.nex tree.%.true.neg.expr.tsv
	/usr/bin/time -o tree.$*.neg.correlation.lrt.lambda.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.neg.expr.tsv \
		--outprefix tree.$*.neg \
		--filter-test lrt \
		--lrt-alt lambda \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run the full LRT phylo signal filter on pos sims
tree.%.pos.correlation.lrt.full.tsv tree.%.pos.correlation.lrt.full.time: tree.%.true.nex tree.%.true.pos.expr.tsv
	/usr/bin/time -o tree.$*.pos.correlation.lrt.full.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.pos.expr.tsv \
		--outprefix tree.$*.pos \
		--filter-test lrt \
		--lrt-alt full \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run the full LRT phylo signal filter on neg sims
tree.%.neg.correlation.lrt.full.tsv tree.%.neg.correlation.lrt.full.time: tree.%.true.nex tree.%.true.neg.expr.tsv
	/usr/bin/time -o tree.$*.neg.correlation.lrt.full.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.neg.expr.tsv \
		--outprefix tree.$*.neg \
		--filter-test lrt \
		--lrt-alt full \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run my implementation of the PATH-based Moran's I autocorrelation phylo signal filter on pos sims
tree.%.pos.correlation.moran.tsv tree.%.pos.correlation.moran.time: tree.%.true.nex tree.%.true.pos.expr.tsv
	/usr/bin/time -o tree.$*.pos.correlation.moran.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.pos.expr.tsv \
		--outprefix tree.$*.pos \
		--filter-test moran \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run my implementation of the PATH-based Moran's I autocorrelation phylo signal filter on neg sims
tree.%.neg.correlation.moran.tsv tree.%.neg.correlation.moran.time: tree.%.true.nex tree.%.true.neg.expr.tsv
	/usr/bin/time -o tree.$*.neg.correlation.moran.time ${GEX_LINEAGE_DIR}/bin/gexLineage \
		--trees tree.$*.true.nex \
		--expr tree.$*.true.neg.expr.tsv \
		--outprefix tree.$*.neg \
		--filter-test moran \
		--filter-only \
		--seed $$(shuf -i 1-1000000000 -n 1)

# Run PATH's implementation of the autocorrelation phylo signal filter on pos sims
tree.%.pos.correlation.path.tsv tree.%.pos.correlation.path.time: tree.%.true.pos.expr.tsv tree.%.true.nex
	singularity exec --bind $(MAIN_DIR):/mnt $(PATH_SIF) /usr/bin/time -o /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.pos.correlation.path.time \
	Rscript /mnt/src/path_gene_phylo_correlation.R \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nex \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.pos.expr.tsv \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.pos.correlation.path.tsv

# Run PATH's implementation of the autocorrelation phylo signal filter on neg sims
tree.%.neg.correlation.path.tsv tree.%.neg.correlation.path.time: tree.%.true.neg.expr.tsv tree.%.true.nex
	singularity exec --bind $(MAIN_DIR):/mnt $(PATH_SIF) /usr/bin/time -o /mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.neg.correlation.path.time \
	Rscript /mnt/src/path_gene_phylo_correlation.R \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.nex \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.true.neg.expr.tsv \
		/mnt/$(REL_PATH)/$(NTAXA)taxa/tree.$*.neg.correlation.path.tsv

# Gather all results to one performance summary file, assuming matched pos and neg results
# Add these files below to the for loop to include all methods
# $(FULLTSVSPOS)
eval.all.performance.txt: $(MORANTSVSPOS) $(MORANTSVSNEG) $(PATHTSVSPOS) $(PATHTSVSNEG) $(LAMBDATSVSPOS) $(LAMBDATSVSNEG) # $(FULLTSVSPOS) $(FULLTSVSNEG)
	{ \
		printf "ntaxa\tmethod\tsim_num\tTP\tFN\tTN\tFP\n"; \
		for f in $(LAMBDATSVSPOS) $(MORANTSVSPOS) $(PATHTSVSPOS); do \
			ntaxa=$(NTAXA); \
			method=$$(basename "$$f" | cut -d"." -f5- | sed 's/\.tsv//' | tr '.' '_'); \
			sim_num=$$(basename "$$f" | cut -d"." -f2); \
			sed -i 's/TRUE/True/g; s/FALSE/False/g' "$$f"; \
			tp=$$(awk -F"\t" 'NR>1 && $$NF=="True" { count++ } END { print count+0 }' "$$f"); \
			fn=$$(awk -F"\t" 'NR>1 && $$NF=="False" { count++ } END { print count+0 }' "$$f"); \
			negf=$$(echo "$$f" | sed 's/\.pos\./\.neg\./'); \
			sed -i 's/TRUE/True/g; s/FALSE/False/g' "$$negf"; \
			tn=$$(awk -F"\t" 'NR>1 && $$NF=="False" { count++ } END { print count+0 }' "$$negf"); \
			fp=$$(awk -F"\t" 'NR>1 && $$NF=="True" { count++ } END { print count+0 }' "$$negf"); \
			echo -e "$$ntaxa\t$$method\t$$sim_num\t$$tp\t$$fn\t$$tn\t$$fp"; \
		done; \
	} > $@

# Gather all runtimes to one summary file, assuming matched pos and neg results
# Add these files below to the for loop to include all methods
# $(FULLTIMESPOS)
eval.all.time.txt: $(MORANTIMESPOS) $(MORANTIMESNEG) $(PATHTIMESPOS) $(PATHTIMESNEG) $(LAMBDATIMESPOS) $(LAMBDATIMESNEG) # $(FULLTIMESPOS) $(FULLTIMESNEG)
	{ \
		printf "ntaxa\tmethod\tsim_num\ttime_sec\n"; \
		for f in $(LAMBDATIMESPOS) $(MORANTIMESPOS) $(PATHTIMESPOS); do \
			ntaxa=$(NTAXA); \
			method=$$(basename "$$f" | cut -d"." -f5- | sed 's/\.time//' | tr '.' '_'); \
			sim_num=$$(basename "$$f" | cut -d"." -f2); \
			time_sec_pos=$$(head -n 1 "$$f" | cut -d" " -f1 | sed 's/user//g'); \
			negf=$$(echo "$$f" | sed 's/\.pos\./\.neg\./'); \
			time_sec_neg=$$(head -n 1 "$$negf" | cut -d" " -f1 | sed 's/user//g'); \
			time_sec=$$(awk "BEGIN {print $$time_sec_pos + $$time_sec_neg}"); \
			echo -e "$$ntaxa\t$$method\t$$sim_num\t$$time_sec"; \
		done; \
	} > $@

eval.all.tree_stats.tsv: $(TREESTATS)
	{ \
		first=1; \
		for f in $(TREESTATS); do \
			ntaxa=$(NTAXA); \
			sim_num=$$(basename "$$f" | cut -d"." -f2); \
			if [ "$$first" = "1" ]; then \
				printf "ntaxa\tsim_num\t"; \
				head -n 1 "$$f"; \
				first=0; \
			fi; \
			tail -n 1 "$$f" | awk -v ntaxa="$$ntaxa" -v sim_num="$$sim_num" 'BEGIN {OFS="\t"} {print ntaxa, sim_num, $$0}'; \
		done; \
	} > $@

clean:
	rm -f tree.* eval.all*

archive-all:
	archive_dir=archive_$(shell date +%Y-%m-%d_%H.%M); \
	mkdir -p $$archive_dir ; \
	mv tree.* $$archive_dir/ ; \
	mv eval.all* $$archive_dir/

