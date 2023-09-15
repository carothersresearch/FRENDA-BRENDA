#!/bin/bash

# Input CSV file and column number (1-based index)
input_file="proteome_exe.csv"
column_number1=1
column_number2=2

# Use 'cut' command to extract the desired column and save it to a variable
column_data1=$(cut -d ',' -f $column_number1 "$input_file")
column_data2=$(cut -d ',' -f $column_number2 "$input_file")

IFS=$'\n' read -d '' -r -a values1 <<< "$column_data1"
IFS=$'\n' read -d '' -r -a values2 <<< "$column_data2"

output_file1="Reaction.csv"
output_file2="SpeciesBaseMechanisms.csv"

# Add header line to the output file
echo "Accession Number,EC,Species,Label,Enzyme,Mechanism,Substrate,Cofactor,Product,Metals,Parameters" > "$output_file1"
echo "Label,Type,StartingConc,Conc,Mechanisms,Parameters" > "$output_file2"

for index in "${!values1[@]}"; do
    value1="${values1[index]}"
    value2="${values2[index]}"
    echo "Processing values: $value1, $value2"

    conversion=$(efetch -db protein -id "$value1" -format gpc -mode xml | xtract -insd Protein EC_number)
    awk -F'\t' -v OFS=',' '{ print $1, $2 }' <<< "$conversion" >> "$output_file1"

    echo ",Enzyme,$value2" >> "$output_file2"
done
