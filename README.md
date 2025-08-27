# BEQ CRISPResso Pipeline
The BEQ (Base Editing Quantifications) CRISPResso Pipeline is a library of functions designed to more efficently perfom CRISPResso functionality on numerous FASTQ files in quick succession. 

CRISPResso_Loop.py recives an amplicon list from the user and runs through large amounts of fastq files directory by directory to eliminate the need for manual file manipulation from a user. 

Requires CRISPResso to be installed on machine.

Upon completion of CRISPRsso_Loop, CRISPResso output is placed in the directory beside coorisponding fastq file(s).

Read_based_quant.py is a python script designed to extract editing data from the CRISPResso output. Specifically a user can select a position within the protospacer corresponding to the intended edit, the user can also select bystanders. The script will give a readout of Allele frequency of the correction both with and without accepted bystander edits and print these to a .csv file.

Upon completion of Read_based_quant, output files are written to a CSV file for the user containing:
- Directory the fastq file was from
- Read-based Lenient Correction Percentage (correction of intended edit along with tolerated bystanders)
- Read-based Strict Correction (correction of intended edit **without** tolerated bystanders)
- Independent Nucleotide Based Editing Quantification Percentage (taken from Quantification_window_nucleotide_percentage_table, allele frequency)
- Reads Aligned
- Reads Total
- Guide Sequence (reverse complement if guide is in reverse orientation from amplicon)
- Guide Sequence relative to amplicon
- Index of Intended Edit
- Index of Permissible Bystanders
- Strings used to search the Alleles_frequency_table in the CRISPResso output

<sub>BEQ CRISPResso Pipeline is distributed under the BSD license.</sub>


## Download ##

BEQ CRISPResso Pipeline can be downloaded from this repository for function use.

Additional information about function implementation can be found in the ```usage``` section

## Directory Format ##
Tree view of working directory before use should display as follows:

<pre><code>
userNamedDirectory/
├── firstDirectoryContainingFastQFiles
│   ├── firstFastqFile.fastq
│   └── optionalAdditionalFastqFile.fastq
├── secondDirectoryContainingFastQFiles
│   └── ...
├── thirdDirectoryContainingFastQFiles
│   └── ...
├── additionalDirectoryContainingFastQFiles
│   └── ...
├── ampliconList.csv
├── CRISPResso_Loop.py
└── Read_based_quant.py
</code></pre>

Tree view of working directory after use will display as follows:

<pre><code>
userNamedDirectory/
├── firstDirectoryContainingFastQFiles
│   ├── firstFastqFile.fastq
│   ├── optionalAdditionalFastqFile.fastq
│   ├── AQLenientCorrection.CSV
│   ├── AQStrictCorrection.CSV
│   ├── CRISPREsso_on_firstFastqFile.html
│   └── CRISPREsso_on_firstFastqFile/
│       ├── 1a.Read_barplot.PDF
│       ├── 1a.Read_barplot.PNG
│       ├── 1b.Alignment_pi_chart.PDF
│       └── ...
├── secondDirectoryContainingFastQFiles
│   └── ...
├── thirdDirectoryContainingFastQFiles
│   └── ...
├── additionalDirectoryContainingFastQFiles
│   └── ...
├── unprocessed_data
│   └── ...
├── ampliconList.csv
├── AQ_Read_Based_Correction
├── CRISPResso_Loop.py
└── Read_based_quant.py
</code></pre>

If more than one fastq file is provided the function will return additional information for additional fastq file.

## SetUp ##
Usage of the pipeline requires the inclusion of Common Amplicon List that contains information about the amplicons in use. The file should be a CSV and should be in the following format:

| Name | Protospacer_sequence | Editor | Guide Orientaion relative to amplicon | Amplicon | note | Tolerated Sequences | Tolerated positions | Intended Edits |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PAH1 | TCA... | ABE | F | CCT... | P281L | TCA... | 16 | 5 |

An example file is included within the respository for reference. If desired, a user will be given the option upon launching the script to create a Common_Amplicon_List. After creation the file it can be placed along side the fastq files awaiting analyzing.

Usage of the pipeline also requires that the fastQ files generated contain the name of their amplicon within the fastq direcotory name. 

## Usage ##

BEQ CRISPResso Pipline can be used in two steps.

### Step 1 ###

First by placing the [CRISPResso_Loop.py](https://raw.githubusercontent.com/Walsh-Quigley/BaseEditingQuantificationsOfCRISPRESSOOutputs/refs/heads/main/CRISPResso_Loop.py)  file in your fastq data directory and then entering 
``` 
python3 CRISPResso_loop.py
```
while in that directory. User will be asked if they'd like to create a Common_amplicon_list file. If needed, the user will have an opportunity to input an amplicon and its relevant information.

After completion the process will exit cleanly.
```
INFO  @ Fri, 27 Jun 2025 01:50:25 (100.0% done):
         Analysis Complete!

INFO  @ Fri, 27 Jun 2025 01:50:25 (100.0% done):

                                        _
                                       '  )
                                       .-'
                                      (____
                                   C)|     \
                                     \     /
                                      \___/
```

CRISPResso output information can be found in each directory containing fastQ files with their respective names.


### Step 2 ###
After the first function exits cleanly enter:
``` 
python3 Read_based_quant.py 
```
User will be prompted for the same information (deliminator and keyword location) and the process will run and exit cleanly. 
```
-----------------
the Directory is: LastDirectoryContaingFastQfile
The Search Term is: XXXXXX
Guide Sequence Variable: XXXXXXXXXXXXXXXXXXXX
Intended Edit/Correction Position: X
Intended Edit/Correction Index (Position -1): X
Guide Orientation: X
Permissible Edit Indicies (Position -1): []
The base at the intended edit index, 8, is an 'X'
Guide is in the XXXXXXX orientation
the corrected guide sequence with intended edit: XXXXXXXXXXXXXXXXXXXX
Found alleles table input file: Alleles_frequency_table_around_sgRNA_XXXXXXXXXXXXXXXXXXXX.txt
New file generated: ../AQLenientCorrection.csv
the lenient correction percentage is: X
Found alleles table input file: Alleles_frequency_table_around_sgRNA_XXXXXXXXXXXXXXXXXXXX.txt
New file generated: ../AQStrictCorrection.csv
the strict correction percentage is: X
The number of reads that align is: XX
The total number of reads: XX
the independent quant sum is: X.XXXXXXXXXXXXX
Current Working Directory: /DATA/LastDirectoryContaingFastQfile
Lenient correction percentage for FastQFileName: X
(base) root@1234abcdef:/DATA#
```
Data returned will be populated with usable information in place of temporary XXXXXXXX.

Comprehensive data compilation can be found in the **AQ_Read_Based_Correction** file in the __userNamedDirectory__ folder

## License ##
BEQ CRISPResso Pipeline is licensed under the BSD 3-clause license.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

The names of its contributors may not be used to endorse or promote products derived from this software without specific prior written permission.


## FAQ ##

## Changelog ##

## Scripts ##

[CRISPResso_Loop.py](https://raw.githubusercontent.com/Walsh-Quigley/BaseEditingQuantificationsOfCRISPRESSOOutputs/refs/heads/main/CRISPResso_Loop.py)