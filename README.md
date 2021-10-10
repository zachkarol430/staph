# STAPH-Gene-Analysis

This code can be used to analyse the staph genome and in the end create a beast phylogenetic tree. A lot of this code is thanks to Rohan so many lines of code, specifically in creating the samples is Rohan's code.

Here are steps to running code(Note this is not command line since many changes may need to be made for each sample)

1. Run create_sample.r and generate desired sample. Do not use too large of a sample as the get_genes function will take a long time. If a large sample is needed then run in batches of 300ish samples and then concatenate samples.
2. Run create_tree.r to create either nexus file which will be used for beast or to create a newick tree.
3. will add extra instruction when clean up files
