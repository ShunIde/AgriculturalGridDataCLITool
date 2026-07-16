# data/widiv_2000SNPs.vcf.gz

Genotype data for 942 lines of the Wisconsin Diversity (WiDiv) maize panel,
2000 SNPs. Vendored unmodified from the [PyBrOpS](https://github.com/rzshrote/pybrops)
example suite:

```
examples/pareto_frontier_visualization/optimal_contribution_selection/widiv_2000SNPs.vcf.gz
```

PyBrOpS is MIT-licensed (Copyright (c) 2019-2023 Robert Shrote). See its
repository for the original license text and the R script used to generate
this file (`widiv_2000SNPs.vcf.gz`'s sibling `*_data_generation.R` in the
same example directory).

`breeding.py` in this repo reads this file with a small stdlib-only parser
rather than PyBrOpS/`cyvcf2` — see the module docstring for why.
