# Variables to add yourself: %TOKEN, %FILENAME, %MODEL
library(staphopia)
library(dplyr)
library(ape)
library(phangorn)
library(phytools)
library(tidytree)
library(ggplot2)
library(ggtree)
library(treeio)


TOKEN = "56166f58c4617fd2c4357794b81cd806d8a3a055" #enter your token
USE_DEV = FALSE
file_name = "final_sample_real" ##change to whatever
setwd("/Users/username/Desktop/Staph")###change to what directory is needed

# useful function
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
# This gets the annotation ids of the core genome
core_index <- read.delim('/Users/zachkarol/iCloud Drive (Archive)/Desktop/Staph/nrd-gene-set.txt', sep = '\t', header=T)
#make this file in directory. Its attached in github

samples <- df$ids[1:20]  ##insert whatever sample you want


# This code block gets the sequence data in nexus format 
t <- proc.time() # keeps track of elapsed time
allg <- get_variant_gene_sequence(as.numeric(samples), annotation_ids = core_index$annotation_id) # gets the core genome sequences for the samples in a non-concatenated format

print(proc.time()-t) # prints elapsed time
allgmr <- subset(allg, sample_id != 'reference') # gets rid of the reference sequence
gallgmr <- allgmr %>% group_by(sample_id) %>% mutate(fullseq  = paste0(sequence, collapse = '')) # concatenates the core genome sequences 
gallgmr <- gallgmr[!duplicated(gallgmr$sample_id),] # I think the previous step has duplicates of samples; this gets rid of any duplicates
nexusdat <- list() 

# this is the data structure to store the sequences
# this loop stores the sequences in the data structure in a way that is understandable by R to write in a nexus format
for(j in 1:length(samples)){
  nexusdat[j] <- strsplit(tolower(gallgmr$fullseq[j]), '')[1]
}
nexusdat<-nexusdat[!is.na(nexusdat)]

names(nexusdat)<- samples

write.table(data.frame(dates=df$date, id=df$ids), file="dates.txt")
write.table(data.frame(id=df$ids,trait=df$mecA),file="trait.txt")#any trait can be added

# adds sample IDs to nexus data
filename <- paste(file_name, '.nex', sep = '')
write.nexus.data(nexusdat, filename) # writes data to filename so you don't have to do this again
print(proc.time()-t)
# Now we want to convert these sequences into phyDat format so that we can make trees.
t <- proc.time()
fbin <- as.DNAbin(nexusdat) # converts the nexus format into DNAbin format (compressed losslessly)
names(fbin) <- names(nexusdat) # names the data 

write.dna(x = fbin, file = paste(file_name, 'dnabin', sep = '')) #writes the DNAbin data to file
fbinphydat <- as.phyDat(fbin) # converts DNAbin data to phyDat data (which is easier to work with)
print(proc.time()-t)
t <- proc.time()
distobj <- dist.dna(fbin, model = "N") # creates distance matrix for samples based on the molecular evolution model %MODEL              # you can try playing around with different options for %MODEL to see how much of a difference it makes
                    fdistmat <- as.matrix(distobj) # turns the "dist" object distobj into a matrix
                    write.csv(x = fdistmat, file = paste(file_name, '_fdistmat.csv', sep = '')) # writes distance matrix to file
                    print(proc.time()-t)
                    # Now finally we make the tree
                    t <- proc.time()
                    njtree <- nj(distobj) # this step performs Neighbor-Joining on the distance matrix to get a tree estimate. This is a commonly-used algorithm but it is important to note that the branch lengths don't mean anything. We use it as a starting point for our tree estimation
                    initpml <- pml(njtree, fbinphydat)
                 # this obtains the initial likelihood for the data under the neighbor joining tree under the default model used by the phangorn package 
                    optimres <- optim.pml(initpml, optNni = TRUE) # this function finds the tree that maximizes the likelihood under the default model used by phangorn, starting with the neighbor-joining tree. It would normally just try to optimize branch lengths, but the optNni option makes it try to optimize topology as well. This will make it much slower, but if it runs in a reasonable time I would strongly recommend it.
                    tr <- optimres$tree
                    write.tree(tr, file = paste(file_name, '.newick', sep = '')) # writes the tree to file in newick format, which is the common format that trees are written in.
                    print(proc.time()-t)
                    # The branch lengths of the resulting tree should be interpreted as time scaled by some kind of population-wide mutation rate.
                    
##depending on purpose newick tree or beast tree can be used
                    

#load beast tree and edit
tree<- read.beast("/Users/zachkarol/Desktop/Staph/new_work.tree")##read in beast file
phy<-multi2di(tree@phylo)
phy$edge.length[phy$edge.length<=0]<-1e-8
tree@data$location.set[tree@data$location.set== "NA"]<- FALSE
tree@phylo<-phy
write.beast(tree, file= "fixed_beast.tree")
##new fixed beast tree which will be used for create_plot script

plotTree(phy,fsize=0.1)#made tip labels small so tree can be seen


















