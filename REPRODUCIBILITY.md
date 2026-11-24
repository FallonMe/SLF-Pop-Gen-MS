# Analysis and reproducibility guide

## Scope

This repository is a research-code archive. The release organizes existing code and metadata and repairs a small number of syntax problems. It does not claim a validated, one-command reproduction of the complete paper.

## Data and execution order

1. Acquire study reads from the two BioProjects linked in the README and the exact study reference assembly. Verify assembly identity before mapping, calling variants, or joining annotations.
2. Review whole-genome workflow notes and adapt the original cluster paths, module names, scheduler resources, and input lists. They include computationally expensive steps; review each section before submitting jobs.
3. Generate the filtered variant and genotype-likelihood inputs required by the downstream analyses. The genotype-likelihood workflow contains separate ANGSD, ngsLD, and PCAngsd stages.
4. Use the relevant `R/` scripts to process outputs, join metadata, and produce figures. They are not numbered execution stages and may reuse in-memory objects or intermediate files from earlier work.
5. Demographic inference requires momi2 and an appropriate site-frequency-spectrum input. Model comparisons and bootstrap summaries require the corresponding fitting outputs.
6. Mitochondrial analyses require their own read/reference inputs, consensus sequences, alignments, and tree/haplotype outputs.

## Environments

- R packages detected in the scripts are listed in `R-packages.txt`. This inventory is not a lockfile; the exact historical installed package versions have not been recovered.
- The Python demographic code imports momi, NumPy, pandas, SciPy, and matplotlib. Consult imports for the individual script/notebook.
- Shell notes reference tools including BWA, samtools, Picard, GATK, vcftools, bcftools, PLINK, ANGSD, ngsLD, PCAngsd, MAFFT, and IQ-TREE. The original notes record some module versions. No portable container is claimed.

## Important execution details

- The extensionless GATK and ANGSD files include prose or multiple languages. Read them as workflow notes; do not pass the entire files to Bash.
- The ANGSD notes contain line-continuation and example-prefix inconsistencies. Review each command against the paper and the intended intermediate filenames before running it. No scientific parameters were silently corrected during packaging.
- Many R files call `setwd()` or contain absolute paths. Adapt them in a separate working copy and send outputs to a new directory. Do not replace the original data or archived results.
- Some scripts depend on objects created interactively. Confirm required inputs and intermediate objects for the section you need.
- Sample identities and ordering must agree across genotype matrices, metadata, coordinates, and environmental matrices. The included validator checks metadata identity consistency; it does not check every genotype file.
- Shared candidate loci and environmental associations should not be read as direct validation of a molecular mechanism.

## Validation performed for this packaging change

- All 21 R scripts copied from `SLF_R_Scripts` pass R syntax parsing after three files received comment-only repairs.
- `momi2_all_model.py` compiles as Python source.
- The two sample metadata tables contain 118 rows each; their IDs and coordinates are checked by `validate_metadata.py`.
- Retained environmental-archive member bytes were preserved while removing Apple/Finder metadata.

The WGS processing, demographic fitting, and complete R analyses have not been rerun locally. Required HPC inputs and exact historical environments are not all present here. No end-to-end reproduction claim is made.
