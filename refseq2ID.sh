#!/bin/bash

# Input CSV file and column number (1-based index)
input_file="proteome_exe.csv"
column_number=2

# Use 'cut' command to extract the desired column and save it to a variable
column_data=$(cut -d ',' -f $column_number "$input_file")

IFS=$'\n' read -d '' -r -a values <<< "$column_data"

output_file="output.csv"

# Add header line to the output file
echo "Accession Number,EC,Species,Name,Substrates,Products,Cofactors,Parameters,Metals,Reaction Type" > "$output_file"


for value in "${values[@]}"; do
    echo "Processing value: $value"
    conversion=$(efetch -db protein -id $value -format gpc -mode xml | xtract -insd Protein EC_number)
    awk -F'\t' -v OFS=',' '{ print $1, $2 }' <<< "$conversion" >> "$output_file"
done
