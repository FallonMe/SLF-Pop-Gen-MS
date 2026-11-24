# R analyses

These scripts were copied from the repository's `SLF_R_Scripts` branch (tip `ffb31022b4c1e69d507052bc861cea2a620e50dc`) so they can be found from the default branch. Existing branches are retained.

| Analysis | Scripts |
| --- | --- |
| Variant summaries | `plot_VCF_features_from_GATKTable.R`, `plot_VCFStats.R` |
| Diversity and homozygosity | `Genome_wide_diversity.R`, `Final_RoH.R` |
| Population structure | `DAPC.R`, `SNMF.R`, `plot_K_admixture.R`, `pcangsd.R` |
| Geographic/environmental distances | `IBD_IBE.R`, `Combine IBD.R` |
| Demographic model summaries | `Momi2_stats.R` |
| Linkage disequilibrium | `LD Decay.R` |
| Nuclear and mitochondrial trees/haplotypes | `SNP_TREE.R`, `Plot_SNPTree.R`, `Plot_MT_tree.R`, `MT_haplotype_analysis.R` |
| Candidate loci and environmental associations | `pcadpt Analysis.R`, `Final_USA_pcadpt.R`, `USA_RDA.R`, `plotRDA.R`, `Shared SNPs het.R` |

All 21 scripts pass R syntax parsing after commenting out Markdown fences and explanatory prose in three scripts. This is a syntax check, not an execution or a statistical reanalysis.

## Running an analysis

These are historical analysis scripts rather than an automated pipeline. Many use original HPC/local paths and intermediate objects. Work through the relevant script and map its inputs to your local data before execution. Do not source every script in sequence or overwrite original study outputs.

`../R-packages.txt` lists packages detected in library()/require() calls; it is not a version lockfile or a complete executable environment. Namespaced calls and command-line tools may require additional dependencies. Computational parameters and statistical methods are retained from the archived scripts.
