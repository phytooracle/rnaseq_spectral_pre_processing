# RNA-Seq & Hyperspectral Pre-processing

## Purpose
This script (i) adds fieldbook data to RNA-Seq and hyperspectral input files, (ii) combines RNA-Seq and hyperspectral data, and (iii) saves the joined RNA-Seq and hyperspectral CSV file to disk for each transcript.

## Arguments
* -r, --rnaseq_csv
    * RNA-Seq input CSV file

* -s, --spectral_csv
    * Hyperspectral input CSV file

* -f, --fieldbook_csv
    * Fieldbook input CSV file

* -o, --out_dir
    * Output directory, default='2019_cotton_rnaseq_spectra'

## Running the code

### Run with input files
You will need to specify the path to each CSV file as such: 

```./rnaseq_spectral_pre_processing.py -r ./data/Cotton_TPM_TOP25.csv -s ./data/2019averagescans.csv -f ./data/Duke_Fieldbook_2019.csv```

### Run without input files
If input files are not specified, the code will automatically read the input files from the CyVerse Data Store using public links:

```./rnaseq_spectral_pre_processing.py```