# Metadata and environmental inputs

Copied from branch `Metadata_EnvData`, tip `5a55fd22a5d8aa3feccbe40c7b3277a9736d0ad7`.

- `Final_SLF_MS_Info.csv`: 118 sample records. Columns: SLF_ID, Country, City, Location, detailed_area, date_col, lat, lon, sex, color, habitat, area.
- `Environmental_PC1s.csv`: 118 sample records with sample information and temp_PC1, precip_PC1, urban_PC1.
- `ENV_Raw_data.zip`: environmental source tables and associated metadata. Apple resource-fork files and Finder metadata were removed from this copy; retained member contents are unchanged.

`lat` and `lon` describe insect collection locations. Field codes and environmental-variable definitions should be interpreted with the published Methods and supplementary materials. The PC scores are the original study values; they are not newly recalculated or interchangeable with later environmental analyses.

Raw sequence data are distributed through NCBI, not GitHub:

- U.S. samples: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1136004
- Shanghai samples: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1258978

Run `python3 validate_metadata.py data` from the repository root to check sample-ID uniqueness, matching sample sets, and coordinate ranges without modifying the data.
