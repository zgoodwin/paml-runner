# paml-runner
This is a set of scripts that determines whether positive selection is acting on a gene by comparing the gene's non-synonymous substitution rate (dN) with its synonymous substitution rate (dS). These rates are calculated by aligning genes from many species and measuring dN and dS. Zhiheng Yang has written a great review on how to interpret these rates \([Open access article here](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2000YangBielawskiTREEv15p496.pdf)\).

[PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) is a suite of tools used to measure whether positive selection is acting on a gene using a statistical test. PAML measures substitution rates (dN and dS) and fits statistical models of positive selection and negative selection to the data. One can then determine which model is more likely using a chi-squared test. It is somewhat tedious to run PAML on more than a few genes at a time. Hence, I wrote this set of scripts to automate the running of PAML on many genes. These scripts were used for an analysis of positive selection in a set of genes required for skin development \([Goodwin et al 2017](http://journal.frontiersin.org/article/10.3389/fgene.2016.00227/full)\). 

For more information on PAML and what it does, you can read its documentation [here](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

## Input file formats

To run this set of scripts, you need a codon alignment (provided in the ```alignments``` folder)

* _Alignments:_ A folder containing alignment files, formatted in the Phylip interleaved format. \([See here.](http://evolution.genetics.washington.edu/phylip/doc/sequence.html)\)
* _Trees:_ Newick format \([See here.](http://evolution.genetics.washington.edu/phylip/newicktree.html)\)
* _Control files:_ The control (or model configuration file) specifies the parameters for each of the PAML models that test for positive selection on a gene. You can create a control file for each PAML model that you would like to run. More information about models and control files can be found in the [PAML documentation](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

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
## parse_paml.py output file format

1. alt_likelihood - The likelihood value of the positive selection model.
2. chi2_stat - The chi-squared statistic (-2\*(Null-alt))
3. dN - The average number of non-synonymous substitutions in the alignment, per alignment site.
4. dN_dS - The average number of non-synoynymous substitutions per synonymous substitution site.
5. dS - The average number of synomnymous substitutions in the alignment, per alignment site. 
6. genes - Gene name
7. null_likelihood - The likelihood value of the negative/neutral selection model.
8. pval - Probability that the null model (negative selection) fits the data. 

## Dependencies

+ Must have a working copy of [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) installed locally
+ scipy
+ numpy

## To-Do's

- [x] Add scripts to repository
- [x] Re-write bash scripts in python
- [x] Re-organize the PAML parser source code
- [x] Update comments
- [x] Update documentation
- [ ] Add plotting scripts
- [ ] Add support for local dN/dS calculation
- [ ] Add support for foreground and background rate calculation
- [ ] Add support for calculating the percentage of sites with a specific dN/dS value
- [ ] Improve error handling for all scripts

