#location of .py file: /Users/aidanq/Desktop/Bash_Scripts/BaseEditingQuantificationsOfCRISPRESSOOutputs/CRISPResso_Loop.py


import os
import subprocess
import glob
import re



# the function to retrieve the GuideSequence AmpliconSequence and GuideOrientation
def retrieveCRISPRessoInputs(search_term):
    guideSequence = ""
    ampliconSequence = ""
    guideOrientation = ""

    with open('./../Common_amplicon_list.csv', 'r') as file:
        for line in file:
            #Split the line by commas and clean up each column
            columns = line.strip().split(',')
            if columns[0].upper() == search_term:
                guideSequence = columns[1].upper().replace("\r", "").replace("-", "")[:20]
                ampliconSequence = columns[4].upper().replace("\r", "").strip()
                guideOrientation = columns[3].upper().replace("\r", "").strip()
                break

    # Print the retrieved variables
    print(f"Guide Sequence Variable: {guideSequence}")
    print(f"Amplicon Sequence Variable: {ampliconSequence}")
    print(f"Guide Orientation: {guideOrientation}")

    return guideSequence, ampliconSequence, guideOrientation

def run_CRISPResso(ampliconSequence, guideSequence, fastq_files):

    if len(fastq_files) == 2:
        subprocess.run([
            'CRISPResso',
            '--fastq_r1', fastq_files[0],
            '--fastq_r2', fastq_files[1],
            '--amplicon_seq', ampliconSequence,
            '--guide_seq', guideSequence,
            '--quantification_window_size', '10',
            '--quantification_window_center', '-10',
            '--base_editor_output'
        ])
    elif len(fastq_files) == 1:
        subprocess.run([
            'CRISPResso',
            '--fastq_r1', fastq_files[0],
            '--amplicon_seq', ampliconSequence,
            '--guide_seq', guideSequence,
            '--quantification_window_size', '10',
            '--quantification_window_center', '-10',
            '--base_editor_output'
        ])
    else:
        print("No fastq files found for this directory")

#find the fastq file names and see if it is single or paired end reads
def gather_fastqs():
        fastq_files = glob.glob('*R1_001.fastq*') + glob.glob('*R2_001.fastq*')
        return fastq_files

def directoryDelimiter():

    # List directories in the current directory
    current_directory = os.getcwd()  # Get the current working cddirectory
    directories = [d for d in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, d))]
    
    # Display the first 3 directories
    print("First 3 directories in the current directory:")
    for i, directory in enumerate(directories[:3], start=1):
        print(f"{i}. {directory}")

    while True:
        #Prompt user for a delimiter
        delimiter = input("Enter the delimiter used in the directories for you data (i.e. did you separate the words/names in the sample sheet by a '-'): ").strip()
        if delimiter.lower() == "\\t":
            delimiter = "\t" #convert the '\t' string to an actual tab character

        if not delimiter:
            print("Delimiter cannot be empty. Please try again")
            continue

        column_input = input("Enter the position in the file name (starting from 1) where the search term is located:").strip()

        if not column_input.isdigit() or int(column_input) < 1:
            print("Invalid column index. Please enter a positive integer.")
            continue

        column_index = int(column_input) - 1
        return delimiter, column_index

# delimiter, column_index = directoryDelimiter()

#Main Loop
for directory in os.listdir():
    if os.path.isdir(directory):
        print("-----------")

        os.chdir(directory)
        print(f"The Directory is: {directory}")

        #gather fastq files
        fastq_files = gather_fastqs()
        print(f"Our fastq file(s) are: {fastq_files}")

        # parts = re.split(delimiter, directory)
        # if len(parts) < 3:
        #     directoryErrorMessage = f"Unexpected directory name format: {directory}"
        #     print(directoryErrorMessage) #print the error
        #     continue #move on to the next directory


        #getting the search term from the directory name using the dash delimiter
        # directoryName = os.path.basename(directory)
        # searchTerm = parts[column_index].upper()
        searchTerm = "AQDB"
        print(f"The Search Term is: {searchTerm}")


        #getting the CRISPResso variables using the search term
        guideSequence, ampliconSequence, guideOrientation = retrieveCRISPRessoInputs(searchTerm)

        #run CRISPResso
        run_CRISPResso(ampliconSequence, guideSequence, fastq_files)

        os.chdir('..')