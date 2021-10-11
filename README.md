# STAPH-Gene-Analysis

This code can be used to analyse the staph genome and in the end create a beast phylogenetic tree. A lot of this code is thanks to Rohan so many lines of code, specifically in creating the samples is Rohan's code.

Here are steps to running code(Note this is not command line since many changes may need to be made for each sample)

1. Run create_sample.r and generate desired sample. Do not use too large of a sample as the get_genes function will take a long time. If a large sample is needed then run in batches of 300ish samples and then concatenate samples.
2. Run create_tree.r to create either nexus file which will be used for beast or to create a newick tree.

3a. If a regular newick tree is created then it will be in your directory and can be loaded into r with read.tree or viewed online with newick viewers
3b. Use nexus file to create beauti xml program. Download Beast 2 and then open Beauti 2. 

4. In Beauti load nexus file and download beast classic(manage packages). Then add tips dates using file "dates" and add discrete trait using file "trait".
5. Fine tune model to liking. I recommend looking at beast website for instructions.
6. Run beast and generate a .trees file and check for convergance in tracer. 
7. Run treeannotator to combine trees and create one target tree.
8. Edit tree in r using create tree file if needed. *Needed if there are NAs or negative branch lengths.
9. Run create_plots to create plots using beast tree. 

Notes
Lots of resources by Liam revell and on beast website. Many plots can be created, I just created a few for simplicity.
