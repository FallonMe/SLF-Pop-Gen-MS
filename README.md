# Population genomics of spotted lanternfly invasion and urbanization

Analysis code and metadata accompanying:

**Meng, F., Snead, A. A., Zhang, A. Y., Munshi-South, J., & Winchell, K. M. (2026). Cities as evolutionary incubators for the global spread of the spotted lanternfly.** *Proceedings of the Royal Society B*, 293, 20252292. [Read the paper](https://doi.org/10.1098/rspb.2025.2292).

This study examines genetic diversity, population structure, demographic history, and candidate genomic responses to environmental variation in the spotted lanternfly (*Lycorma delicatula*) across native Shanghai and introduced U.S. populations.

## Start here

- **Understand the analyses:** read the paper and the [analysis guide](REPRODUCIBILITY.md).
- **Find the R code:** browse [R/](R/README.md).
- **Find sample/environmental inputs:** browse [data/](data/README.md).
- **Check the metadata:** run the command below.
- **Cite the work:** use [CITATION.cff](CITATION.cff).

```bash
python3 validate_metadata.py data
```

This check uses Python's standard library and does not run genomic analyses or change input files.

## Repository map

| Location | Contents |
| --- | --- |
| `R/` | Population structure, diversity, demographic summaries, phylogenetic visualization, environmental association, and candidate-locus analyses |
| `data/` | Sample metadata, original study environmental PC scores, and environmental source tables |
| `Momi2_models.ipynb`, `momi2_all_model.py` | Demographic model exploration and fitting with momi2 |
| `SLF_GATK_HPC_Pipeline_Final` | Historical whole-genome processing notes, including shell and Python sections |
| `ANGSD_pipeline` | Historical ANGSD/ngsLD/PCAngsd workflow notes |
| `Mt Dna Phylogeography Pipeline hpc` | Mitochondrial processing, consensus, alignment, and phylogenetic workflow |
| `REPRODUCIBILITY.md` | Dependencies, workflow order, expected inputs, and execution limits |
| `CHANGELOG.md` | Documentation and packaging changes |

## Computational workflow

1. Process whole-genome reads, align to the study reference, and call/filter variants; examine the GATK workflow notes.
2. Use genotype-likelihood analyses where appropriate; examine the ANGSD workflow notes.
3. Analyze diversity, population structure, geographic/environmental associations, and candidate loci with the R scripts.
4. Fit and compare demographic models using the momi2 code and summarize their outputs in R.
5. Process mitochondrial sequences and visualize phylogenetic/haplotype results.

These are historical research scripts, not a single executable workflow. Some guides mix prose and languages; many analyses depend on study-specific paths or intermediate files. Follow the [reproducibility notes](REPRODUCIBILITY.md) before running them. Syntax checks do not establish that the full study has been reproduced.

## Data access

The study includes 98 U.S. and 20 Shanghai whole-genome samples. Raw sequencing reads are available through NCBI BioProjects [PRJNA1136004](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1136004) and [PRJNA1258978](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1258978), respectively. Reference-genome context: [Snead, Meng et al. 2025](https://doi.org/10.1038/s41597-025-04854-8). Use the exact reference assembly and filtering settings described in the study; do not silently substitute a different assembly.

## Provenance and versions

The `SLF_R_Scripts` and `Metadata_EnvData` branches are retained as historical records. Their relevant files are also available on the default branch. Packaging changes and limited syntax repairs are recorded in [CHANGELOG.md](CHANGELOG.md); no statistical models or study results were refitted for this release.

## Credit and contact

Study authors: Fang Meng, Anthony A. Snead, Aria Yang Zhang, Jason Munshi-South, and Kristin M. Winchell. External tools and packages retain their authorship and licenses; cite the original methods described in the paper.

Questions: [Fang (Fallon) Meng](mailto:fm2371@nyu.edu). [Research portfolio](https://fallonme.github.io/).
