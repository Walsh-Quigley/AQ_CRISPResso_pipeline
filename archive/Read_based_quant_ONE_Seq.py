#location of .py file: /Users/aidanq/Desktop/Bash_Scripts/BaseEditingQuantificationsOfCRISPRESSOOutputs/Read_based_quant.py

import os
import re
from itertools import combinations
import csv
import glob
from datetime import datetime
from pathlib import Path
import shutil

# Global variables
analysis_result_file = "../AQ_Read_Based_Correction.csv"    

#function to get all our variables from the common amplicons list
def gimmieDat(searchTerm):
    guideSequence = ""
    guideOrientation = ""

    with open('./../Common_amplicon_list.csv', 'r') as file:
        match_found = False
        for line in file:
            #split the line by commas and clean up each column
            columns = line.strip().split(',')
            if columns[0].upper() == searchTerm:
                guideSequence = columns[1].upper().replace("\r", "").replace("-", "")[:20]
                guideOrientation = columns[3].upper().replace("\r", "").strip()
                match_found = True
                break

        if not match_found:
            raise ValueError(f"Search term '{searchTerm}' not found in file.")


    # Print the retrieved variables
    print(f"Guide Sequence Variable: {guideSequence}")
    print(f"Guide Orientation: {guideOrientation}")

    #checking to see if the guide is in the forward orientation and calling reverse complement function if not
    if guideOrientation == 'F':
        print("Guide is in the forward orientation")
    else:
        print("the guide is in the reverse orientation, generating reverse complement...")
        guideSequence = reverseComplement(guideSequence)

    return guideSequence, guideOrientation



def reverseComplement(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Replace each base with its complement and reverse the resulting sequence
    return ''.join(complement[base] for base in reversed(dna))

def generate_toleratedSequences(sequence, indices, replacement_letter):
    toleratedSequences = []

    # Iterate over all possible numbers of positions to replace (1 to len(indices))
    for r in range(1, len(indices) + 1):
        # Generate all combinations of r indices
        for combo in combinations(indices, r):
            # Create a new sequence for this combination of indices
            new_sequence = list(sequence)  # Convert to list for mutability
            for index in combo:
                new_sequence[index] = replacement_letter
            toleratedSequences.append(''.join(new_sequence))  # Convert back to string

    return toleratedSequences

def CRISPRessoDirectoryHelperFunction():
    # Get a list of all items in the current directory
    current_dir = os.getcwd()
    subdirs = [d for d in os.listdir(current_dir) if os.path.isdir(d) and d.startswith("CRISPResso_")]

    if not subdirs:
        print("No directory starting with 'CRISPResso_' found.")
        return False

    # Assuming you want to change to the first matching directory
    target_dir = subdirs[0]
    os.chdir(target_dir)
    return True

def sum_last_column(output_file):
    total = 0
    with open(output_file, 'r') as outfile:
        reader = csv.reader(outfile)
        next(reader)  # Skip the header row
        for row in reader:
            try:
                # Convert the last column to a float and add it to the total
                total += float(row[-1])  # Assuming the last column is numeric
            except ValueError:
                # Handle case where the last column is not numeric
                print(f"Skipping row with non-numeric value in the last column: {row}")
                continue

    return total

def append_header_and_timestamp(csv_file):
    # Check if the file exists
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, 'a', newline='') as result_file:
        writer = csv.writer(result_file)

        # Add header if the file doesn't exist
        if not file_exists:
            writer.writerow(['Directory',
                             'Percent Edited (any "A")', 
                             'Reads Aligned', 
                             'Reads Total', 
                             'Guide Sequence ("-" removed, truncated to 20bp)',
                             'Guide Sequence relative to amplicon',
                             'Index of "A"s ',
                             'Strings used to search the Alleles_frequency_table in the CRISPResso output'])  # Write the header

        # Get the current date and time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Append the last run line with timestamp
        writer.writerow([f"the below analysis was run @ {current_time}"])

def allelesTableFilter(directory, arrayOfStrings, output_file):
    if type(arrayOfStrings) == 'str':
        print(f"Incorrect Data Structure passed: {type(arrayOfStrings)}")

    correctionSum = -1
    CRISPRessoDirectoryHelperFunction()

    input_files = glob.glob("Alleles_frequency_table*.txt")

    if not input_files:
        print("No file starting with 'Alleles_frequency_table' found.")
        correctionSumError = "No file starting with 'Alleles_frequency_table' found."
        os.chdir("..")
        errorPrintStatement(directory, correctionSumError)
        return correctionSum
    
    else:
        input_file = input_files[0]  # Take the first matching file
        print(f"Found alleles table input file: {input_file}")

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')  # Tab-delimited reader
        writer = csv.writer(outfile)  # Default is comma-delimited for CSV

        #move the header to the output file
        header = next(reader)
        writer.writerow(header)

        # Process each row
        for row in reader:
            row_first_column = row[0].strip().upper()
            if any(s.strip().upper() in row_first_column for s in arrayOfStrings):  # Check if any string in the array matches
                writer.writerow(row)  # Write matching row to the output file

    correctionSum = sum_last_column(output_file)

    print(f"New file generated: {output_file}")

    os.chdir("..")

    return correctionSum

def extractReadCounts():
    readsAligned = -1
    readsTotal = -1

    CRISPRessoDirectoryHelperFunction()

    CRISPRessoMappingStatsFile = "CRISPResso_mapping_statistics.txt"


    file_path = Path(CRISPRessoMappingStatsFile)

    if file_path.is_file():
        with file_path.open("r") as file:
            header = file.readline()
            data_line = file.readline().strip()
            col = data_line.split("\t")
            readsAligned = int(col[2])
            readsTotal = int(col[0])
    else:
        print(f"Error: File '{CRISPRessoMappingStatsFile}' does not exist.")
        os.chdir("..")
        readsAligned = "The 'CRISPResso_mapping_statistic.txt' file does not exist"
        readsTotal = "The 'CRISPResso_mapping_statistic.txt' file does not exist"
        return readsAligned, readsTotal

    print(f"The number of reads that align is: {readsAligned}")
    print(f"The total number of reads: {readsTotal}")

    os.chdir("..")

    return readsAligned, readsTotal

def generateQuantificationOutput(directory, 
                                 lenient, 
                                 readsAligned, 
                                 readsTotal, 
                                 guideSequence,
                                 guideOrientation,
                                 permissibleEdits,
                                 toleratedSequences):
        
        directoryTruncated = directory.split("_L001")[0]
        analysis_result_file = "../AQ_Read_Based_Correction.csv"

        print(f"Lenient correction percentage for {directoryTruncated}: {lenient}")
        with open(analysis_result_file, 'a', newline='') as result_file:
            writer = csv.writer(result_file)
            writer.writerow([directoryTruncated, lenient, readsAligned, readsTotal, guideSequence, guideOrientation, permissibleEdits, toleratedSequences])

def independentQuant(correctionLocationIndex, guideOrientation):
    independentQuantSum = -1

    CRISPRessoDirectoryHelperFunction()
    nucleotidePercentageFile = "Quantification_window_nucleotide_percentage_table.txt"

    if guideOrientation == "R":
        row_index = 2
    else:
        row_index = 3

    col_index = int(correctionLocationIndex) + 1

    file_path = Path(nucleotidePercentageFile)

    if file_path.is_file():
        with open(nucleotidePercentageFile, "r") as file:
            lines = file.readlines()

            target_row = lines[row_index].strip()
            col = target_row.split("\t")
            independentQuantSum = col[col_index]
    else:
        print(f"Error: File '{nucleotidePercentageFile}' does not exist")
        os.chdir("..")
        independentQuantSum = "Error: the file 'Quantification_window_nucleotide_percentage_table.txt' does not exist"
        print(independentQuantSum)
        errorPrintStatement()
    
    independentQuantSum = float(independentQuantSum) * 100

    os.chdir("..")
    return independentQuantSum

def directoryDelimiter():

    # List directories in the current directory
    current_directory = os.getcwd()  # Get the current working cddirectory
    directories = [d for d in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, d))]
    
    # Display the first 3 directories
    print("First 3 directories in the current directory:")
    for i, directory in enumerate(directories[:3], start=1):
        print(f"{i}. {directory}")

    while True:
        # Prompt user for delimiters
        delimiter_input = input(
            "Enter the delimiter(s) used in the directories for your data (e.g., '-', '_'). "
            "For multiple delimiters, enter them without spaces (e.g., '-_'): "
        ).strip()
        
        if not delimiter_input:
            print("Delimiter cannot be empty. Please try again.")
            continue

        # Convert input into a regex pattern to match any of the provided delimiters
        delimiter_pattern = f"[{re.escape(delimiter_input)}]"

        column_input = input("Enter the position in the file name (starting from 1) where the search term is located:").strip()

        if not column_input.isdigit() or int(column_input) < 1:
            print("Invalid column index. Please enter a positive integer.")
            continue

        column_index = int(column_input) - 1
        return delimiter_pattern, column_index

def errorPrintStatement(directory, directoryErrorMessage):
     
    with open(analysis_result_file, 'a', newline='') as result_file:
            writer = csv.writer(result_file)
            writer.writerow([directory, directoryErrorMessage])

    print("there appears to be an error in processing this directory")
    moveFile = False
    currentDirectory = os.getcwd()
    print(currentDirectory)
    userMoveFileResponse = input("Would you like to generate a directory called 'unprocessed_data' and move this file to it (y or n)").lower()
    
    if userMoveFileResponse == "y":
        moveFile = True

    if moveFile == True:
        print("Moving this file to the unprocessed_data directory")
        parent_directory = os.path.dirname(currentDirectory)
        unprocessed_data_directory = os.path.join(parent_directory, "unprocessed_data")

        if not os.path.exists(unprocessed_data_directory):
            os.makedirs(unprocessed_data_directory)

        target_path = os.path.join(unprocessed_data_directory, os.path.basename(currentDirectory))
        shutil.move(currentDirectory, target_path)
        print(os.getcwd())
        os.chdir("..")

    print(os.getcwd())

def identifyPermissibleEdits(guideSequence):

        permissibleEditsIndices = [i for i, letter in enumerate(guideSequence) if letter == "A"]
        return permissibleEditsIndices








delimiter, column_index = directoryDelimiter()

append_header_and_timestamp("AQ_Read_Based_Correction.csv")

#rewriting the main loop
for directory in os.listdir():
    #checks if the current item (directory) is a directory and not a file
    #if it returns true then the next block of code is allowed to run 
    if directory == "unprocessed_data":
        print("Skipping this directory: {directory}")
        continue

    if os.path.isdir(directory):
        print("-----------------")
        print(f"the Directory is: {directory}")

        #changes the working directory to the current directory meaning 
        #subsequent operations will be relative to this subdirectory
        os.chdir(directory)

        parts = re.split(delimiter, directory)
        if len(parts) < 3:
            directoryErrorMessage = f"Unexpected directory name format: {directory}"
            print(directoryErrorMessage) #print the error
            errorPrintStatement(directory, directoryErrorMessage) #log the error
            continue #move on to the next directory


        #getting the search term from the directory name using the dash delimiter
        directoryName = os.path.basename(directory)
        searchTerm = parts[column_index].upper()
        print(f"The Search Term is: {searchTerm}")

        #passing that search term to the function that will extract all the info from 
        #common amplicons list
        try:
            guideSequence, guideOrientation = gimmieDat(searchTerm)
        except Exception as e:
            print(f"Error in common amplicons list: {searchTerm}: {e}")
            errorPrintStatement(directory, e)
            os.chdir("..")
            continue

        permissibleEditsIndicies = identifyPermissibleEdits(guideSequence)

        #generating the sequences that we will filter for
        toleratedSequences = generate_toleratedSequences(guideSequence, permissibleEditsIndicies, "G")

        #search the alleles freq table for our tolerated sequences
        lenientCorrectionFile = "../AQLenientCorrection.csv"
        lenientCorrectionSum = allelesTableFilter(directory, toleratedSequences, lenientCorrectionFile)
        if lenientCorrectionSum == -1:
            os.chdir("..")
            continue
        print(f"the lenient correction percentage is: {lenientCorrectionSum}")

        #search alleles freq table for the strict edit
        # correctedSequenceArray = [correctedSequence]
        # strictCorrectionFile = "../AQStrictCorrection.csv"
        # strictCorrectionSum = allelesTableFilter(directory, correctedSequenceArray, strictCorrectionFile)
        # print(f"the strict correction percentage is: {strictCorrectionSum}")

        #grabbing the read counts
        readsAligned, readsTotal = extractReadCounts()

        # #grabbing independent editing quantification data from quantification_window_nucleotide_percentage_table.txt
        # independentQuantSum = independentQuant(correctionLocationIndex, guideOrientation)

        # print(f"the independent quant sum is: {independentQuantSum}")


        print("Current Working Directory:", os.getcwd())

        #generate the table containing the quantification output
        generateQuantificationOutput(directory, lenientCorrectionSum, readsAligned, readsTotal, guideSequence, guideOrientation, permissibleEditsIndicies, toleratedSequences)

        os.chdir("..")

