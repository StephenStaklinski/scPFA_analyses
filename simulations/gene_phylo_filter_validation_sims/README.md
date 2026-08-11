To run all simulations, set `THREADS` in the top-level Makefile and then use:

```bash
make run-all
```

To generate tree statistics for each simulation size, use:

```bash
make run-treestats
```

To collect the performance, runtime, and tree-statistic results and generate all
plots, activate an environment containing NumPy, pandas, Matplotlib, seaborn, and
SciPy, then use:

```bash
make collect-results
```
