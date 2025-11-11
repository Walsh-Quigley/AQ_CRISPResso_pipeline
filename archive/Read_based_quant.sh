# Function to extract guide sequence, amplicon sequence, guide orientation, and edits
extract_guide_info() {
    echo "Current Directory: $DIR"
    local searchTerm="$1"

    # Convert searchTerm to uppercase for case-insensitive matching
    searchTerm=$(echo "$searchTerm" | tr '[:lower:]' '[:upper:]')
    
    # Extract guideSeq, ampSeqVar, guideOrientation, and edit indices
    guideSeq=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $2}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | tr -d '-' | xargs | cut -c1-20)
    ampSeqVar=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $5}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)
    guideOrientation=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $4}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)
   
    intendedEditIndex=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $9}' ./../Common_amplicon_list.csv | tr -d '\r' | xargs)
    permissibleEditIndex=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $8}' ./../Common_amplicon_list.csv | tr -d '\r' | xargs)
    
    permissibleEditArray=($permissibleEditIndex)

    echo "Extracted Intended Edit Index (before conversion): $intendedEditIndex"

    # Validate the intendedEditIndex before attempting conversion
    if [[ -z "$intendedEditIndex" || ! "$intendedEditIndex" =~ ^[0-9]+$ ]]; then
        echo "Error: intendedEditIndex is empty or not a valid integer."
        return 1
    fi

    # Convert intendedEditIndex to an integer if valid
    intendedEditIndex=$(printf "%d" "$intendedEditIndex")
    echo "Intended Edit Index (after conversion): $intendedEditIndex"
}

# Function to check if the character at intendedEditIndex in guideSeq is "A"
check_intended_edit() {

    local guideSeq="$1"
    local intendedEditIndex="$2"

        # Check if intendedEditIndex is a valid integer
    if [[ -z "$intendedEditIndex" || ! "$intendedEditIndex" =~ ^[0-9]+$ ]]; then
        echo "Error: intendedEditIndex is empty or not a valid integer."
        return 1
    fi

    local index=$(($intendedEditIndex - 1))  # Convert to 0-based index

    local charAtIndex=${guideSeq:$index:1}

    if [[ "$charAtIndex" == "A" ]]; then
        echo "The character at position $intendedEditIndex in guideSeq is A."
    else
        echo "The character at position $intendedEditIndex in guideSeq is not A, it's $charAtIndex."
    fi
}

# Function to generate all combinations of permissible edits
generate_combinations() {
    local baseSeq="$1"
    local intendedEditIndex="$2"
    local permissibleEditArray=("${!3}")

    # Create an array to hold all possible sequences with combinations of permissible edits
    local result=()

    # Start with the intended edit (replace 'A' at intendedEditIndex with 'G')
    local index=$((intendedEditIndex - 1))  # Convert to 0-based index
    local firstEdit="${baseSeq:0:$index}G${baseSeq:$((index + 1))}"
    result+=("$firstEdit")

    # Generate all combinations of permissible edits
    local numPermissible=${#permissibleEditArray[@]}
    local numCombinations=$((2 ** numPermissible))

    for (( i=1; i<numCombinations; i++ )); do
        modifiedSeq="$firstEdit"  

        for (( j=0; j<numPermissible; j++ )); do
            if (( (i >> j) & 1 )); then
                permissibleIndex=$((permissibleEditArray[j] - 1))
                charAtPermissibleIndex=${modifiedSeq:$permissibleIndex:1}

                # Only replace 'A' with 'G' at permissible positions
                if [[ "$charAtPermissibleIndex" == "A" ]]; then
                    modifiedSeq="${modifiedSeq:0:$permissibleIndex}G${modifiedSeq:$((permissibleIndex + 1))}"
                fi
            fi
        done

        result+=("$modifiedSeq")
    done

    echo "${result[@]}"


}

# Function to generate reverse complement of a DNA sequence
reverse_complement() {
    local seq="$1"
    # Reverse the sequence and replace the nucleotides with their complements
    echo "$seq" | rev | tr 'ATCG' 'TAGC'
}


# Main loop: editing quantification loop
for DIR in */; do

    echo "--------------"

    # Move into the directory
    cd "$DIR";

    # Get the directory name
    directoryName=$(basename "$DIR")

    # Turn the name into a search term for searching the spreadsheet
    searchTerm=$(echo "$directoryName" | awk -F'[-_]' '{print $2}')

    # Print the search term for visibility
    echo "The search term is: $searchTerm"
    
    # Extract relevant variables using the function
    extract_guide_info "$searchTerm"

    # Print the array elements
    echo "The silent or tolerated bystanders are at positions: ${permissibleEditArray[@]}"

    # Call the check_intended_edit function, which handles the check and print logic
    check_intended_edit "$guideSeq" "$intendedEditIndex"

    # Generate all combinations of permissible edits
    permissibleEditStrings=()

    permissibleEditStrings=($(generate_combinations "$guideSeq" "$intendedEditIndex" permissibleEditArray[@] | tail -n 1))
    echo "Permissible edit combinations generated: ${#permissibleEditStrings[@]} sequences."

    # Handle guideOrientation being "R" for reverse complement
    if [[ "$guideOrientation" == "R" ]]; then
        echo "Guide orientation is R. Generating reverse complements..."
        for (( i=0; i<${#permissibleEditStrings[@]}; i++ )); do
            reverseComp=$(reverse_complement "${permissibleEditStrings[$i]}")
            permissibleEditStrings[$i]="$reverseComp"
            echo "Reverse complement of ${permissibleEditStrings[$i]}: $reverseComp"
        done
    else
    echo "The guide is in the same orientation as the amplicon."
    fi

    # Print the final array of modified sequences
    echo "Final search terms with permissible edits: ${permissibleEditStrings[@]}"

    # Navigate to the CRISPResso folder and generate the CSV file from the table
    cd CRISPR*
    readBasedAlignmentTxt=$(ls Alleles_frequency_table_around_sgRNA_*.txt 2>/dev/null)
    readBasedAlignmentTable="Alleles_frequency_table_around_sgRNA.csv"
    sed 's/\t/,/g' "$readBasedAlignmentTxt" > "$readBasedAlignmentTable"
    lenientCorrection="corrected_lenient.csv"

    # Extract total reads and reads aligned from CRISPResso_mapping_statistics.txt
    totalReads=$(awk 'NR==2 {print $1}' ./CRISPResso_mapping_statistics.txt)
    readsAligned=$(awk 'NR==2 {print $3}' ./CRISPResso_mapping_statistics.txt)

    # Print total reads and reads aligned for visibility
    echo "Total Reads: $totalReads"
    echo "Aligned Reads: $readsAligned"


    # Track if any matching rows are found
    foundMatch=false

    for string in "${permissibleEditStrings[@]}"; do
        # Use a temporary variable to store the output and check if awk matches any line
        awkResult=$(awk -v search="$string" -F',' '$1 ~ search {print $0}' "$readBasedAlignmentTable")

        # Only append if awkResult is not empty
        if [[ -n "$awkResult" ]]; then
            echo "$awkResult" >> "$lenientCorrection"
            foundMatch=true
        fi
    done


    # Check if no matches were found; if so, print 0
    if ! $foundMatch; then
        lenientCorrectionSum=0
        echo "LINE 146: NO READS ALIGNED"
    else
        # Sum the last column of the output file (excluding the header)
        lenientCorrectionSum=$(awk -F',' 'NR>1 { total += $NF } END { print total }' "$lenientCorrection")
        echo "Lenient correction sum: $lenientCorrectionSum"
    fi


    cd ..

    # Combine searchTerm with extracted values and write to CSV
    final="$directoryName,$lenientCorrectionSum,$totalReads,$readsAligned"
    echo "$final" >> ./../Correction_Read_Based.csv

    # Move back out of the directory to the main directory
    cd ..;

done

#add a header to our table
echo -e "directoryName,lenientCorrectionPercentage,totalReads,readsAligned" | cat - Correction_Read_Based.csv > temp && mv temp Correction_Read_Based.csv
