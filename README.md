# Analysis scripts for Ohigashi _et al._ (2025) _Soil Biol. Biochem._
These scripts are provided to reproduce the results and figures presented in Ohigashi _et al._ (2025) _Soil Biology & Biochemistry_. DOI: [10.1016/j.soilbio.2025.109895](https://doi.org/10.1016/j.soilbio.2025.109895)

The preprint is also avalable here: [![bioRxiv](https://img.shields.io/badge/Preprint-bioRxiv-red.svg)](https://doi.org/10.1101/2025.02.17.638669)

---

## Important Notes
### 1. Microbial count data
The raw microbial count data (`rarefied_ASV_table_16S.txt` and `rarefied_ASV_table_ITS.txt`) are not included in this repository, although references to them are made in the scripts.  

To generate these files, please follow the steps below:

1. Download the single-end FASTQ files from NCBI SRA:  
   - **PRJNA699079** (for soils collected in Kenya)  
   - **PRJNA699085** (for soils collected in Malawi)  
2. A correspondence table between NCBI accession numbers (`Run`) and FASTQ file names is provided in `Data/Accession_num_list.xlsx`.  
3. Perform quality filtering, denoising, and creation of ASV tables using DADA2. Use the following scripts for these steps:  
   - `01_1_DADA2_16S.R`  
   - `01_2_DADA2_ITS.R`  
   Running these scripts will generate `rarefied_ASV_table_16S.txt` and `rarefied_ASV_table_ITS.txt`.

### 2. Microbial function analysis
To execute `02_2_PICRUSt2.sh`, you need to install PICRUSt2 on your system beforehand. Detailed installation instructions can be found [here](https://github.com/picrust/picrust2).

