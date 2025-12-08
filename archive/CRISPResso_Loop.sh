# Function to retrieve guideSeqVar, ampSeqVar, and guideOrientation
get_variables() {
    local searchTerm=$1
    guideSeqVar=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $2}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | tr -d '-' | xargs | cut -c1-20)
    ampSeqVar=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $5}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)
    guideOrientation=$(awk -F',' -v searchTerm="$searchTerm" 'toupper($1) == searchTerm {print $4}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)

    echo "Guide Sequence Variable: $guideSeqVar"
    echo "Amplicon Sequence Variable: $ampSeqVar"
    echo "Guide Orientation: $guideOrientation"
}

# Function to run CRISPResso
run_crispresso() {
    local ampSeq=$1
    local guideSeq=$2
    local fastqFiles=(*.fastq*)

    if [[ ${#fastqFiles[@]} -eq 2 ]]; then
        CRISPResso \
        --fastq_r1 *R1_001.fastq* \
        --fastq_r2 *R2_001.fastq* \
        --amplicon_seq "$ampSeq" \
        --guide_seq "$guideSeq" \
        --quantification_window_size 10 \
        --quantification_window_center -10 \
        --base_editor_output
    else
        CRISPResso \
        --fastq_r1 *R1_001.fastq* \
        --amplicon_seq "$ampSeq" \
        --guide_seq "$guideSeq" \
        --quantification_window_size 10 \
        --quantification_window_center -10 \
        --base_editor_output
    fi
}

# Main Loop
for DIR in */; do

    echo "----------"

    # Move into the directory
    cd "$DIR";

    # Get the directory name
    directoryName=$(basename "$DIR")

    # Turn the name into a search term for searching the spreadsheet, converting to uppercase
    searchTerm=$(echo "$directoryName" | awk -F'[-_]' '{print $2}' | tr '[:lower:]' '[:upper:]')

    # Print the search term for visibility in the console
    echo "Search Term: $searchTerm"

    # Get guideSeqVar, ampSeqVar, and guideOrientation
    get_variables "$searchTerm"

    # Run CRISPResso with the retrieved variables
    run_crispresso "$ampSeqVar" "$guideSeqVar"

    # Move back out of the directory to the main directory
    cd ..;

done
