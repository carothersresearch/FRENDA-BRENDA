#!/bin/bash

if [ -z "$1" ]; then
  filtoption="unfiltered"  # Set a default value if the third argument is not provided
else
  filtoption="$1"  # Use the provided value for the filtoption
fi

# Call the Python script with the input file as argument
python parser.py output.csv brenda_download.txt "$filtoption"
