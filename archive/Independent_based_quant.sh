
#editing quanification loop
for DIR in */; do

    echo "------------"

    #move into the directory
    cd "$DIR";

    #get the directory name
    directoryName=$(basename "$DIR")

    #turn the name into a search term for searching the spreadsheet
    searchTerm=$(echo "$directoryName" | awk -F'[-_]' '{print $2}')


    #Print the search term so we can have some visibility in the console
    echo "searchterm: $searchTerm"
    
    #grab relevant variables
    guideSeqVar=$(awk -F',' -v searchTerm="$searchTerm" '$1 == searchTerm {print $2}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | tr -d '-' | xargs | cut -c1-20 )
    ampSeqVar=$(awk -F',' -v searchTerm="$searchTerm" '$1 == searchTerm {print $5}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)
    guideOrientation=$(awk -F',' -v searchTerm="$searchTerm" '$1 == searchTerm {print $4}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)
    intendedEditIndex=$(awk -F',' -v searchTerm="$searchTerm" '$1 == searchTerm {print $9}' ./../Common_amplicon_list.csv | tr '[:lower:]' '[:upper:]' | tr -d '\r' | xargs)

    #Print these variables so we have visibility in the console
    echo "guideSeqVar: $guideSeqVar"
    echo "ampSeqVar: $ampSeqVar"
    echo "guideOrientation: $guideOrientation"
    
    #create an empty array to store the positions of our relevant character
    positions=()

    #----------------------This should be an if else
    #starts with the assumption of searching for 'A's in the forward orientation
    targetChar="A"
    sedCommand='4!d'

    echo "line 78: $sedCommand"

    #changes this assumption if the guide is in the reverse orientation relative to the amplicon
    if [[ "$guideOrientation" == "R" ]]; then
        targetChar="T"
        sedCommand='3!d'
    fi
    #----------------------
    #else: you broke it


    echo "line 85: $sedCommand"

    #this grabs the columns of the quantification table and turns it into a string because this is how CRISPResso's output works
    editingWindow=$(head -n 1 ./*/quantification_window_nucleotide_percentage_table.txt | tr -d '\t')

    #we then loop through the editing window seeing which characters match our target character, 'A' for ABEs in the forward orientation
    for ((i=1; i<=${#editingWindow}; i++)); do
        char="${editingWindow:i-1:1}"        
        
        if [[ "$char" == "$targetChar" ]]; then
            positions+=($i)
        fi
    done

    #prints in the terminal the positions that are recognized as matching
    echo "Positions of '${targetChar}' after for loop: ${positions[@]}"

    #turns this into text to later add to a column in a spreadsheet for more visibility to users of this code
    positionsText=$(echo "Positions of '${targetChar}': ${positions[@]}")

    echo "Intended Edit location: $intendedEditIndex"

    #generates another array to store values of each position previously identified
    extracted_values=()


    # Loop through positions array, extract values, and store them
    for pos in "${positions[@]}"; do
        value=$(awk -v col="$pos" '{print $(col+1)}' ./*/quantification_window_nucleotide_percentage_table.txt | sed "$sedCommand")
        extracted_values+=("$value")
    done

    echo "extracted values: ${extracted_values[*]}" 

    #print out the contents of CRISPResso_mapping_statistics.txt
    totalReads=$(awk 'NR==2 {print $1}'  ./*/CRISPResso_mapping_statistics.txt)
    readsAligned=$(awk 'NR==2 {print $3}'  ./*/CRISPResso_mapping_statistics.txt)

    echo "Total Reads: $totalReads"
    echo "Aligned Reads: $readsAligned"



    # Combine searchTerm with extracted values and write to CSV
    final="$directoryName,$searchTerm,$totalReads,$readsAligned,$guideOrientation,$targetChar,$guideSeqVar,$editingWindow,$positionsText,$intendedEditIndex,$(IFS=,; echo "${extracted_values[*]}")"
    echo "$final" >> ./../Independent_Editing_Frequency.csv


    #move back out of the directory to the main directory
    cd ..;

done

#add a header to our table
echo -e "directoryName,AmpliconName,totalReads,readsAligned,guideOrientation,targetCharacter,guideSequence,editingWindow,positions,intendedEditIndex" | cat - Independent_Editing_Frequency.csv > temp && mv temp Independent_Editing_Frequency.csv
