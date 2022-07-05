# MeDIP

## Description
This repository for the MeDIP-seq data analysis workflow in R.
This workflow built on the [QSEA package](http://bioconductor.org/packages/release/bioc/html/qsea.html) and can analyze DNA methylation at gene levels. 

## Instructions
- The **prepare_samples** folfer contain the bash script to prepare the data (from raw .fastq file) before doing the QSEA analysis in R
- The **QSEA_Rcode_in_server** and **QSEA_Rcode_server_2021Jan18** are the R script to run QSEA analysis.
- The **data** folder stores all the needed data information and output after the QSEA analysis.
- The **R_code_forACM**, **R_code_forFrontier**, and **R_code_forFrontier_RNAseq** are the R script to run the downstream analysis for the publications in ACM and Frontiers in Bioscience - Landmark journal, respectively.
- The **outcome_ACM** and **outcome_Frontiers** folders store all the outcome for the publications in ACM and Frontiers in Bioscience - Landmark journal, respectively.

Other files: 
- **R_functions** - store all the functions that are used in other R scripts
- **anotation_NN_2020Nov16** - script to prepare the annotation files

## Related publications 
Nhan Nguyen, Matthias Lienhard, Ralf Herwig, Jos Kleinjans, and Danyel Jennen. “Epirubicin alters DNA methylation profiles related to cardiotoxicity”. Frontiers in Bioscience - Landmark (2022, in press).

Nhan Nguyen, Matthias Lienhard, Ralf Herwig, Jos Kleinjans, and Danyel Jennen. “A bioinformatics workflow to detect genes with DNA methylation alterations: a case study of analyzing MeDIP-seq data in cardiac microtissue exposed to epirubicin”. International Conference Proceedings by ACM (2022). [https://doi.org/10.1145/3510427.3510437](https://doi.org/10.1145/3510427.3510437)

## Contact:
LinkedIn:	[nhannguyen](https://www.linkedin.com/in/nhannguyen1412) | ORCID: [0000-0001-8720-1195](https://orcid.org/0000-0001-8720-1195)
