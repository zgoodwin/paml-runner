!#/bin/bash

# Set up the analysis directories
python setup_paml.py ./control_files/master_M8.ctl ./alignments/ ./trees/14_species_great_ape.nwk M8
python setup_paml.py ./control_files/master_M7.ctl ./alignments/ ./trees/14_species_great_ape.nwk M7

# Run PAML's codeml algorithm
python run_paml.py M8 /Users/zanegoodwin/Desktop/paml4.8/bin/codeml
python run_paml.py M7 /Users/zanegoodwin/Desktop/paml4.8/bin/codeml

# Parse the results and calculate p-values
python parse_paml.py H1/ H0/ m7_m8_results.txt