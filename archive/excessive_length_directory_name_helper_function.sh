#!/bin/bash
#file location: /Users/aidanq/Desktop/Bash_Scripts/excessive_length_directory_name_helper_function.sh

# Loop through directories starting with ds.
for DIR in ds.*; do
    # Check if the directory contains a .fastq.gz file with L001 in the name
    fastqFile=$(find "$DIR" -type f -name '*L001*.fastq.gz' -print -quit)

    if [[ -n "$fastqFile" ]]; then
        # Extract the part of the file name before L001
        fileName=$(basename "$fastqFile")
        truncatedFileName=${fileName%%L001*}

        # Extract the directory name
        dirName=$(basename "$DIR")

        # Combine the truncated file name with the directory name
        newDirName="${truncatedFileName}-${dirName}"

        # Rename the directory
        echo "Renaming directory: $DIR -> $newDirName"
        mv "$DIR" "$newDirName"
    else
        echo "No matching .fastq.gz file found in $DIR"
    fi
done
