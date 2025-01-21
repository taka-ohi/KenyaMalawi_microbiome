# R and Shell scripts for Ohigashi et al. (202X) _Journal Name_
These scripts are to reproduce results and figures in Ohigashi et al. (202X) _Journal Name_. Citation information is as follows: _coming soon_ 
<br>
<br>

## Important Notes
### 1. Microbial count data
- The raw microbial count data ("rarefied_ASV_table_16S.txt" and "rarefied_ASV_table_ITS.txt") are not shared in this repository although they appear in the scripts.
- Before running the scripts, you need to download our single-end fastq files from NCBI SRA XXXXXX. There is a correspondence table for the NCBI accession numbers and file names in "Data/Accession_num_list.csv".
- Then, to do quality filtering, denoizing, and making ASV table, you need to run DADA2 using codes in "01_1_DADA2_16S.R" and "01_2_DADA2_ITS.R", which produces "rarefied_ASV_table_16S.txt" and "rarefied_ASV_table_ITS.txt".

### 2. Microbial function analysis
- To run "02_2_PICRUSt2.sh", you need to install PICRUSt2 on your PC first. The instruction is [here](https://github.com/picrust/picrust2/wiki/Installation).

