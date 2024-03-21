## The CURED Pipeline: Classification Using Restriction Enzyme Diagnostics

---

## **Introduction**
While whole genome sequencing (WGS) is the gold standard for surveying such pathogens, access to resources and the associated costs make applying WGS for larger numbers of isolates prohibitive in resource-rich settings, and impossible in resource-scarce settings. Consequently, the need for bioinformatics tools that can leverage sequencing data to quickly find diagnostic sequence biomarkers that can identify the molecular epidemiology of clones and outbreaks, is growing. To close this gap, we present the **Classification Using Restriction Enzymes Diagnostics (CURED)** pipeline, which uses a limited amount of local sequencing data and data from public databases to output all of the information required for a restriction-enzyme based diagnostic test to identify a specific clone of interest. This novel pipeline promises to be a catalyst for equitable disease surveillance worldwide, transforming the way developing public health risks are confronted and contained.


---
## **Installation**

### Dependencies

- mlst (tested with version 2.23.0)
- ncbi-datasets (tested with version 15.28.0)
- unitig-caller (tested with version 1.3.0)
- samtools (tested with version 1.18)
- blast (tested with version 2.15.0)
- bwa (tested with version 0.7.17)
- biopython (tested with version 1.82)

#### Installing Dependencies

```conda create --name cured```

```conda activate cured```

```conda install -c conda-forge -c bioconda -c defaults mlst==2.23.0```

```conda install -c conda-forge ncbi-datasets-cli==15.28.0```

```conda install -c bioconda unitig-caller```

```conda install -c bioconda blast==2.15.0```

```conda install -c bioconda samtools```

The CURED pipeline was tested with Linux and macOS.

---
## **CURED Toolbox**
(1) **CURED_Main.py**: The purpose of this script is to identify any sequences that are unique to the designated case group. The sensitivity and specificity thresholds can be defined by the user. The default is 100% sensitivity and 100% specificity.  
(2) **CURED_FindREs.py**: The purpose of this script is to identify any unique restriction enzyme sites in any of the unique k-mers. The specificity for finding these sites can be set by the user. By default, specificity is set to 100%, meaning that by finding the same enzyme site in just one control precludes the restriction site from being unique.  

---
### Running CURED_Main.py:

This script has multiple entry points:

##### (1) Using a curated dataset with case/control CSV file.
Formatting of the CSV file should look as follows:

GCA_XXXXXXX.X.fna,control  
GCA_XXXXXXX.X.fna,case

By default, the program expects the extension .fna. Use the --extension flag to indicate a different extension.

Running CURED_Main.py:  
``` python3 CURED_Main.py --case_control case_control_input.csv --genomes_folder /path/to/genomes/```

 Running CURED_Main.py with a curated dataset but with altered sensitivity and specificity thresholds:  
 ``` python3 CURED_Main.py --sensitivity 95 --specificity 95 --case_control case_control_input.csv --genomes_folder /path/to/genomes/```

 Running CURED_Main.py with a curated dataset but with k-mer size of 31 (the default is 20-mer):  
 ``` python3 CURED_Main.py --kmer_length 31 --case_control case_control_input.csv --genomes_folder /path/to/genomes/```  


##### (2) Selecting a species and designating a specific sequence type as the case group.
 ``` python3 CURED_Main.py --species Staphylococcus aureus --sequence_type 300```

 By default, this will download all *Staphylococcus aureus* assembled genomes available in NCBI. Use the --database option to choose to download genomes from just RefSeq or GenBank.

 ```python3 CURED_Main.py --species Staphylococcus aureus --sequence_type 300 --database RefSeq```

 You can also see how many genomes are to be downloaded by ncbi-datasets before running the above command by using the ```--summary``` option.

 ```python3 CURED_Main.py --species Staphylococcus aureus --sequence_type 300 --database RefSeq --summary```  

 #### (3) Providing a list of accession numbers to be used as the case group.
Formatting of the file:  
GCA_000013125  
GCA_000439755  
GCA_000439775  

 ```python3 CURED_Main.py --species Pseudomonas putida --case_accession_list case_accessions.txt --database GenBank```  

 #### (4) Providing a list of accession numbers to be used as cases and controls.
 Formatting of the file:  
 GCA_000013125,case  
 GCA_000439755,case  
 GCA_000439775,control  
 GCA_000478425,control  

 ```python3 CURED_Main.py --case_control_file case_control_accessions.txt --use-datasets```  

 #### (5) Providing a folder of local genomes to be used as cases and download a species to be used as controls.

 ```python3 CURED_Main.py --case_genomes --genomes_folder genomes/ --species Anaplasma phagocytophilum --database GenBank```  

#### (6) Providing a list of known k-mers that can be queried against a set of genomes.  
If you have a set of k-mers that you are interested in, you can leverage the ```--use_simple``` option to search a set of genomes for the presence of these specific k-mers. This option is different than the previous options in that it skips the first iteration of unitig-caller (which is run in call mode) and goes directly into running unitig-caller in simple mode. This option can be used with any of the previous entry point options.

Formatting of the file:  
kmer   
AAAAAAAAAAAAAAAA  
AATTTTTTTAAAGGGCCCCC

*Note the kmer header of this file.*

Using simple mode with a curated dataset:  

```python3 CURED_Main.py --use_simple --kmer_list kmers_to_query.txt --case_control_file case_control.txt --genomes_folder /path/to/genomes/ --extension fa```    

Using simple mode with a species and ST of interest:   

```python3 CURED_Main.py --use_simple --kmer_list kmers_to_query.txt --species Staphylococcus aureus --sequence_type 105 --database RefSeq```

```
usage: CURED_Main.py [-h] [--number_of_cases NUMBER_OF_CASES] [--number_of_controls NUMBER_OF_CONTROLS]
                     [--threads THREADS] [--case_control_file CASE_CONTROL_FILE] [--sensitivity SENSITIVITY]
                     [--specificity SPECIFICITY] [--kmer_length KMER_LENGTH] [--species SPECIES [SPECIES ...]]
                     [--sequence_type SEQUENCE_TYPE] [--case_accession_list CASE_ACCESSION_LIST] [--extension EXTENSION]
                     [--database DATABASE] [--summary] [--quiet] [--genomes_folder GENOMES_FOLDER] [--case_genomes]
                     [--use_datasets] [--use_simple] [--kmer_list KMER_LIST]

This script is part of the CURED pipeline. This script is used for finding unique clonal biomarkers in your cases.

options:
  -h, --help            show this help message and exit
  --number_of_cases NUMBER_OF_CASES, -cases NUMBER_OF_CASES
                        Add in the number of cases to be used in each iteration.
  --number_of_controls NUMBER_OF_CONTROLS, -controls NUMBER_OF_CONTROLS
                        Add in the number of controls to be used in each iteration.
  --threads THREADS, -T THREADS
                        Number of threads to be used for unitig-caller. Default = 1.
  --case_control_file CASE_CONTROL_FILE
                        Csv file of genomes to be used with case or control designation.
  --sensitivity SENSITIVITY
                        Specifiy sensitivity. Default = 100.
  --specificity SPECIFICITY
                        Specify specificity. Default = 100.
  --kmer_length KMER_LENGTH, -K KMER_LENGTH
                        Specify minimum length of k-mer to search for. Default is 20.
  --species SPECIES [SPECIES ...]
                        Genus or species of interest
  --sequence_type SEQUENCE_TYPE, -st SEQUENCE_TYPE
                        sequence type of interest
  --case_accession_list CASE_ACCESSION_LIST
                        List of case accessions.
  --extension EXTENSION, -x EXTENSION
                        extension of assembly inputs. Ignore if using --species/--sequence_type options. Default = fna
  --database DATABASE, -db DATABASE
                        Choose to download genomes from RefSeq, GenBank, or both. Default = both.
  --summary             Check to see how many genomes will be downloaded if you use the --species option. Use this
                        option with --database and --species options.
  --quiet, -q           No screen output. Default = OFF
  --genomes_folder GENOMES_FOLDER
                        Path to genomes.
  --case_genomes        Option if you have local sequencing data to serve as the cases. Use --species to download
                        control genomes.
  --use_datasets        Option to provide CURED with a case and control file of ncbi accessions and downloaded them. Use
                        with --case_control_file
  --use_simple          Option to run unitig-caller simple mode. This is useful for when you already have a list of
                        k-mers that you want to query against a set of genomes.
  --kmer_list KMER_LIST
                        List of k-mers to be used as query in running unitig-caller simple mode
```
#### Output from running CURED_Main.py
 File | Description
 ---- | -----------
Unique_Kmers_Report.txt | Report of unqiue k-mers identified within set sensitivity and specificity. Reports the number of controls and/or cases that each k-mer is found in.
tmp*/ | Intermediate raw data files produced and used by unitig-caller
UniqueKmers.txt | List of unique k-mers. This can be used as input for CURED_FindREs.py
*.log | Log report

### Running CURED_FindREs.py
By default, after identifying the k-mer in a case genome, the program adds 20 bases to either end of the sequence. This is to make it discernible on a gel once amplified. The number of bases added to each end can be adjusted using the --add_bases flag.

The *UniqueKmers.txt* file from running CURED_Main.py can be used directly as input in this script.

#### Examples of running CURED_FindREs.py:

```python3 CURED_FindREs.py --case_control_file case_control_input.csv UniqueKmers.txt /path/to/genomes/```

Change the specificity of each restriction enzyme site  
```python3 CURED_FindREs.py --case_control_file case_control_input.csv --specificity 50 UniqueKmers.txt /path/to/genomes/```

To find the support for each restriction enzyme site found in the k-mer, set --specificity 0  
```python3 CURED_FindREs.py --case_control_file case_control_input.csv --specificity 0 UniqueKmers.txt /path/to/genomes/```

The default mode for identifying unique restriction enzyme sites is based upon presence/absence. For instance, the position in which the restriction enzyme is found in within the control sequence is not taken into account. If the same restriction enzyme is found in the control as in the case, it's not considered unique to the case. However, the mode can be changed using the ```--compare_coordinates``` flag. When used, this flag will consider the position of the restriction enzyme site in the control. If the same restriction enzyme in the case is found in the flanking region of the control, then it is considered to be unique to the case.  

```python3 CURED_FindREs.py --compare_coordinates UniqueKmers.txt /path/to/genomes/```

Use ```--compare_coordinates``` option with ```--specificity 0``` to find the full support for each identified restriction enzyme site using comparison mode.

```python3 CURED_FindREs.py --compare_coordinates --specificity 0 UniqueKmers.txt /path/to/genomes/```


```
usage: CURED_FindREs.py [-h] [--case_control_file CASE_CONTROL_FILE]
                        [--specificity SPECIFICITY] [--added_bases ADDED_BASES]
                        [--min_coverage MIN_COVERAGE] [--extension EXTENSION]
                        [--pcr_product_upstream PCR_PRODUCT_UPSTREAM]
                        [--pcr_product_downstream PCR_PRODUCT_DOWNSTREAM]
                        [--compare_coordinates]
                        kmers genomes_folder

This script is a part of the CURED pipeline. This script is used to find restriction enzyme
sites in the identified k-mers.

positional arguments:
  kmers                 List of kmers to be searched.
  genomes_folder        Path to genomes.

optional arguments:
  -h, --help            show this help message and exit
  --case_control_file CASE_CONTROL_FILE
                        Csv file of cases and controls.
  --specificity SPECIFICITY, -S SPECIFICITY
                        Specificity for finding RE sites in controls. Default = 100.
  --added_bases ADDED_BASES
                        Number of added bases on either end of the sequence. Default = 20.
  --min_coverage MIN_COVERAGE
                        Minimum coverage threshold for sequence to be considered found in
                        controls. Default = 90.
  --extension EXTENSION, -x EXTENSION
                        Extension of input assemblies. Default = fna
  --pcr_product_upstream PCR_PRODUCT_UPSTREAM, -UP PCR_PRODUCT_UPSTREAM
                        Number of bases to include upstream of identified k-mer in outputted
                        PCR product.
  --pcr_product_downstream PCR_PRODUCT_DOWNSTREAM, -DOWN PCR_PRODUCT_DOWNSTREAM
                        Number of bases to include downstream of identified k-mer in
                        outputted PCR product.
  --compare_coordinates
                        Mode to compare RE by position to determine uniqueness. Default is to
                        determine uniqueness based on presence/absence.
```
#### Output from running CURED_FindREs.py
 File | Description
 ---- | -----------
CURED_UniqueEnzymes.tsv | Report of unique restriction enzyme sites in the k-mers.
CURED_FindREs_controls.txt | If a restriction enzyme site that is found in the case genomes is identified in any of the control genomes, the names of the controls are reported with the corresponding k-mer and restriction enzyme.
CURED_FindREs_summary.txt | Summary report detailing the k-mer, case genome used, number of control genomes used in the search and the number of control genomes excluded due to no alignment found, or low coverage.
CURED_UniqueEnzymes_PCR_Products.tsv | If a restriction enzyme is identified as unique to the case, the PCR product for the corresponding k-mer is outputted.
*.log | Log report
