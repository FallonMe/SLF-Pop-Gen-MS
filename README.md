This repository contains all scripts, pipelines, and metadata used in our study, “Cities as evolutionary incubators for the global spread of the spotted lanternfly”. This study investigates the role of urbanization in shaping the genetic diversity, population structure, demographic history, and local adaptation of the invasive spotted lanternfly (Lycorma delicatula).

We integrate whole genome resequencing, mitochondrial phylogeography, and eco-evolutionary modeling to dissect how cities facilitate rapid evolution and global biological invasions.

1. Momi2_models.ipynb/
   
Contains all demographic inference analyses performed with momi2, including:
Construction of the Site Frequency Spectrum (SFS)
Eight alternative demographic scenarios
Visualization of demographic history and effective population sizes

3. SLF_R_scripts/
   
This script will contain all R scripts used for:
Genetic diversity statistics (ROH,heterozygosity etc.)
Population structure (PCA, Admixture, sNMF)
Isolation by Distance / Environment (IBD, IBE)
Genotype-environment association (RDA)
Selection scan integration and GO term enrichment (pcadapt etc.)

5. SLF_Analysis_Pipeline_HPC/
   
High-performance computing (HPC) shell scripts for preprocessing and analysis of whole-genome data using tools available on NYU’s Greene HPC cluster.
These include pipelines for:
Quality control: FastQC
Read alignment: BWA-MEM
Mark duplicates: Picard
Variant calling and filtering: GATK (HaplotypeCaller, GenotypeGVCFs, hard-filtering)
Mitochondrial genome mapping: bwa mem, samtools, bcftools, and custom consensus scripts
Post-processing: vcftools, bcftools, plink, and others
Refer to SLF_GATK_HPC_Pipeline_Final/ and ANGSD_pipeline/ for complete example job scripts.

7. ANGSD_pipeline/
   
Contains scripts and parameter files used with ANGSD for:
Estimating genome-wide diversity
Population structure and genetic distance.

9. Mt_DNA_Phylogeography_Pipeline/
    
Complete pipeline for mitogenomic haplotype analysis:
Consensus mitogenome generation
Multiple sequence alignment using MAFFT
Maximum likelihood phylogenetic reconstruction with IQ-TREE
Includes scripts for integrating our samples with published mitogenomes (Du et al. 2021)

Metadata Branch
Metadata for 118 SLF samples and Environmental Data:
Contains metadata files for the 118 spotted lanternfly samples used in the study, along with all associated environmental data. 

For any questions or issues related to this repository, please contact fm2371@nyu.edu.
