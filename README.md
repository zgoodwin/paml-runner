# paml-runner
Determines whether positive selection is acting on one or more genes by comparing their non-synonymous substitution rate (dN) with their synonymous substitution rate (dS). These rates are calculated by aligning genes from many species and measuring dN and dS.

[PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) is a suite of tools used to measure whether positive selection is acting on genes and proteins. Unfortunately, it is difficult to run PAML on more than a few genes at a time. Hence, I wrote this set of scripts to automate the running of PAML on many genes. These scripts were used for an analysis of positive selection in a set of genes required for skin development \([Goodwin et al 2017](http://journal.frontiersin.org/article/10.3389/fgene.2016.00227/full)\). 

## Example

When running PAML, one normally fits codon substitution data for multiple species to two statistical models. One model assumes that positive selection is acting on a gene (dN/dS > 1), while another assumes that negative, or purifying selection is acting on the gene (dN/dS <= 1). Briefly, here's how to do this with paml-runner:

```bash
# Set up the analysis directories to fit codon data to a postitive selection model (M8) and a negative selection model (M7)
python setup_paml.py ./control_files/master_M8.ctl ./alignments/ ./trees/14_species_great_ape.nwk M8
python setup_paml.py ./control_files/master_M7.ctl ./alignments/ ./trees/14_species_great_ape.nwk M7

# Run PAML's codeml algorithm on the same data for each of the two models
python run_paml.py M8 /Users/zanegoodwin/Desktop/paml4.8/bin/codeml
python run_paml.py M7 /Users/zanegoodwin/Desktop/paml4.8/bin/codeml

# Parse the results and calculate p-values
python parse_paml.py M8/ M7/ m8_m7_results.txt
```

## Dependencies

+ Must have a working copy of [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) installed locally
+ scipy
+ numpy

## To-Do's

- [x] Add scripts to repository
- [x] Re-write bash scripts in python
- [x] Re-organize the PAML parser source code
- [ ] Add plotting scripts
- [ ] Add support for local dN/dS calculation
- [ ] Add support for foreground and background rate calculation
- [ ] Add support for calculating the percentage of sites with a specific dN/dS value
- [ ] Improve error handling for all scripts
- [x] Update comments
- [x] Update documentation
