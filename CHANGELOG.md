# Changelog

## Public-release preparation — 2026-11-24

- Added a project README with publication, workflow, data-access, and citation information.
- Copied 21 R scripts from `SLF_R_Scripts` into `R/`, retaining the source branch.
- In `IBD_IBE.R`, `Momi2_stats.R`, and `SNP_TREE.R`, commented out Markdown fences and explanatory prose that prevented R parsing. Computational expressions and statistical settings were not changed.
- Copied the two sample metadata tables from `Metadata_EnvData` into `data/`, retaining the source branch.
- Repackaged `ENV_Raw_data.zip` without Apple resource-fork or Finder files; retained scientific data bytes are unchanged.
- Added dependency inventory, metadata validation, citation metadata, and execution limitations.
- Preserved the original main-branch analysis files and notebook outputs.

This release does not refit analyses or alter the paper's scientific conclusions. The more recent scripts found locally were reviewed separately and were not substituted for the GitHub archive.
