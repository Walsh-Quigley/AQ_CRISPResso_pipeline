# AQ CRISPResso Pipeline 

## Introduction
This pipeline allows the automation of CRISPResso analysis on fastq files with minimum amount of user input as well as creating a summary of the relevant information created by the CRISPResso function. 

## Prerequisites
This pipeline is created for Mac or Linux use. If working on a windows machine it is recommended that a user install a Linux distribution such as Ubuntu that can be used through WSL. The following installation guide will be written for implementation within a Conda environment. The list of software needed is listed here. A detailed installation guide will follow.
- Conda/Miniconda
- Python 3.13 or later
- CRISPResso 2 
- Pandas (Python Package)

## Installation
### Step 1: Install Conda
Conda is a package and environment manager that is recommended for use if the user is working with multiple pipelines that need multiple different environments. This installation will follow the path of Miniconda which only includes Conda and Python. The download link for Miniconda can be found [here](https://www.anaconda.com/download/success), download the version fit for your machine. It is of note that if installing on Windows, you will also need to download a WSL interface as well (Ubuntu is recommended, which can be found in the Microsoft Store).

### Step 2: Creating and activating a Conda Environment
This step will be a guide to creating a Conda environment. Open the terminal on your machine. For Mac users this will be Terminal, for Windows users this will be Ubuntu. Next we will make a Conda environment that has Python installed using this command, the name of the environment can be chosen by the user:
  ```bash
  conda create -n [Name Of Environment] python=3.13
  ```
We can then activate the Conda environment to begin analysis:
  ```bash
  conda activate [Name Of Environment]
  ```
When working with this pipeline within Conda we will need to use this command to activate the environment. 

### Step 3: Installing CRISPResso2 and Python Packages
We will now install CRISPResso2 using this command:
  ```bash
  conda install -c bioconda crispresso2
  ```

To check that the package has been installed correctly we can use the command:

  ```bash
  CRISPResso --version
  ```
Which will return the user's CRISPResso version. We will also install Pandas for our quantification summary using this command:
  ```bash
  pip install pandas
  ```

### Step 4: Downloading the script itself
Click on the Code button near the top of the page to download the ZIP file containing the main functions along with their dependent scripts. This ZIP can be moved and extracted to the user's desired location. This location will be where we perform the data analysis. Your fastq files should also be in this directory, within a folder named "fastqs". A directory tree will be pictured below.

### Step 5: Navigating to the project directory
Returning to the user's terminal, we will now navigate to the project directory. On MacOS the command will look similar to:

```bash
cd path/to/YourProjectDirectoryName
```
On Ubuntu the path will look similar to:
```bash
/mnt/c/Users/YourName/YourProjectDirectoryName
```
Your project directory should resemble this:
```
YourProjectDirectoryName/
├── CRISPResso_Loop.py          # Main script for Pipeline 1
├── Quantification_loop.py       # Main script for Pipeline 2
├── amplicon_list.csv            # Configuration file 
├── scripts/                     # Helper functions folder
│   ├── verify_amplicon_list.py
│   ├── identify_amplicon.py
│   ├── run_CRISPResso.py
│   └── ... (other helper scripts)
└── fastqs/                      # Your sequencing data goes here
    ├── Sample1/
    │   ├── reads_R1.fastq.gz
    │   └── reads_R2.fastq.gz
    ├── Sample2/
    │   └── reads.fastq
    └── ...
```
## Configuring amplicon_list.csv
The final file required before analysis is the amplicon_list CSV. This CSV will contain the following:
- name
- protospacer_or_PEG
- editor
- guide_orientation_relative_to_amplicon
- amplicon
- note
- tolerated_edits
- intended edits 

### Example amplicon_list.csv

| name | protospacer_or_PEG | editor | guide_orientation_relative_to_amplicon | amplicon | note | tolerated_edits | intended_edit |
|------|---------------------|--------|----------------------------------------|----------|------|-----------------|---------------|
| PAH1 | TCA... | ABE | F | CCt... | | 16 | 5 |
| ATL1 | GTC... | ABE | F | TTT... | | | ONE-seq |


A sample amplicon_list.csv is available within the zip file downloaded. This file can be altered by the user to fit their needs. It is of note that the file does not have to be named "amplicon_list" exactly but the file name MUST contain the phrase "amplicon_list".
Examples of acceptable names:
- amplicon_list
- 2025_amplicon_list
- amplicon_list_01012025

Examples of unacceptable names:
- amplicons
- amplicon_2025_list
- amplicon list

Relevant notes: If performing ONEseq analysis, instead of writing a number in the intended edit column, the user can simply write "ONE-seq".

# Pipeline 1: CRISPResso Loop
## Overview
The first main script is the CRISPResso_Loop.py. This script allows for CRISPResso analysis to be performed on all files within your fastq directory in quick succession with minimal user input.

## 1: Directory Structure, Amplicon List, and Naming Convention
 Each fastq subdirectory within the "fastqs" directory must contain the name of the of the amplicon you intend to use as your CRISPResso input. For example, if you want the PAH1 amplicon to be used for a fastq directory, the directory must contain PAH1 in it, as such:
```
YourProjectDirectoryName/
└── fastqs/                      # Your sequencing data goes here
    ├── Sample1_PAH1_directory/
    │   └── reads_R1.fastq.gz
```
The number of .gz files within the fastq subdirectory should never exceed 2, though the pipeline is equipped to deal with both single and double end reads.

## 2: Executing CRISPResso_Loop.py
After verifying the contents of your amplicon list and directory, as well as making sure you are in the working directory and (Conda environment if applicable) for your analysis, we are now ready to run the script itself. The script can be run using this command:
```bash
python CRISPResso_Loop.py
```
After the command is executed, the terminal will display the analysis progress as it moves through the fastq subdirectories. Run times will vary based on the number of fastq files in the sequence, generally run times for each subdirectory are between 15 and 60 seconds. Once analysis has completed, the terminal will display a brief summary indicating if any subdirectories were skipped or if any errors were encountered. 

## 3: Understanding Output
### Logfiles
Upon executing these scripts, log files are generated describing the run in more depth and indicating what errors, if any, were encountered. If an error is encountered based on directory format or an invalid amplicon list, the likely cause will be written to the log file. In the case of a failed run or unexpected values after analysis, it is suggested that the user to check the log file as many errors can be remedied with minimal effort from the user.

### CRISPResso Output
The analysis from the CRISPResso function itself can be found within the corresponding subdirectory that it refers to. Within each subdirectory will be the original unmodified .gz files, the CRISPResso output HTML, as well as a subdirectory containing all of the files and images that make up the output HTML. Format for a correctly analysed fastq sub directory should take this form:
```
YourProjectDirectoryName/
└── fastqs/                      # Your sequencing data goes here
    ├── Sample1_PAH1_directory/
        ├── reads_R1.fastq.gz
        ├── CRISPResso_on_reads_R1.html
        └── CRISPResso_on_reads_R1/
            ├── 1a.Read_barplot.PDF
            ├── 1a.Read_barplot.PNG
            ├── 1b.Alignment_pie_chart.PDF
            └── ...
```
All information procured from CRISPResso analysis can be found here. In addition to this, the subsequent quantification loop uses these subdirectories to create the summary for the user. For that reason CRISPResso_Loop.py must always be executed prior to Quantification_Loop.py. 

# Pipeline 2: Quantification Loop
## 1: Overview
The quantification loop generates a extensive summary of the information outputted by CRISPResso. This script provides data point for both regular and ONE-seq analysis. Because the directory structure will match that of the CRISPResso_Loop.py section, we will move to execution and output directly. 

## 2: Executing Quantification_loop.py
Similarly to CRISPResso_Loop.py, once the user verifies that their directory format is correct and they are in their desired environment, the script can be executed using this command:
```bash
python Quantification_loop.py
```
The terminal will once again display the progress as the data is parsed and compiled. The run time will be in ratio to the first scripts run time but should take slightly less time. The user will be asked if they would like to generate data formatted for a PRISM table which can be confirmed or declined. Finally a brief summary will be displayed indicating any issue or errors encountered during the run.

## 3: Understanding Output
### Logfiles
Log files will once again be created for the user after a run, indicating any issues and offering potential solutions to common errors. Log files will be stored in a directory named "logs" within the users working directory.

### Quantification loop output:
#### ONEseq output
All information will be contained within CSV files placed in the fastqs directory. 
If a user intends to produce ONE-seq analysis from the script the output table will contain the following information:
| Column | Header | Description |
|--------|--------|-------------|
| A | Sample | Fastq directory identifier |
| B | Reads aligned | Number of sequencing reads that align with associated amplicon |
| C | Reads total | Total number of sequencing reads in the fastq files |
| D | Percentage of reads with A>G in first 10 base pairs | Percentage of reads showing ANY A -> G edits in the first 10 base pairs of the given protospacer sequence |
| E | Percentage of reads with A>G in protospacer | Percentage of reads showing A -> G edits at any point within the given protospacer |
| F | Guide sequence | |
| G | A>G 10 base pairs search sequence | All the search sequences checked that contain an A -> G edit in the first 10 base pairs |
| H | A>G any search sequence | All possible A -> G edit combinations in the entire protospacer |
#### Non-ONEseq output
If a user intends to produce analysis containing intended and tolerated edit positions, the the output table will contain the following information:
| Column | Header | Description |
|--------|--------|-------------|
| A | Sample | Fastq directory identifier |
| B | Reads aligned | Number of sequencing reads that align with associated amplicon |
| C | Reads total | Total number of sequencing reads in the fastq files |
| D | Correction with bystanders | Percentage of reads where correction exists at intended position, allowing for intended bystanders |
| E | Correction without bystanders | Percentage of reads where correction exists at intended position, with no other edits |
| F | Independent correction | Percentage of reads where an A is edited to a G anywhere within the protospacer |
| G | Independent correction MINUS correction with bystanders | |
| H | Correction with bystanders MINUS correction without bystanders | Percentage of reads with tolerated bystander edits |
| I | Target locus | Guide sequence used for a given sample |
| J | Perfect correction | Ideal edited sequence with only the intended correction |
| K | Corrected locus with bystander | All edited sequences searched |

If a user's amplicon table contains both ONEseq and non-ONEseq analysis then two distinct tables will be produced.


