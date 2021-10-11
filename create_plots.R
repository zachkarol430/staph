library("runner")
library(plotly)


tree<- read.beast("/Users/zachkarol/Desktop/Staph/it_now_works.tree")
phy<-multi2di(tree@phylo)

plotTree(phy, node.numbers=T)#super messy


#may be useful if tips get messed up and needed to recalculate genes
genes<-get_genes(as.numeric(unlist(tree@phylo$tip.label)),exclude_sequence = T)
allgenesn <- subset(genes, name != 'none')
genetable <- table(allgenesn$sample_id, allgenesn$name)
genetable[which(genetable > 1)] <- 1
df<-as.data.frame.matrix(genetable) 
df$date<-as.Date(sample_id_list_meta$collection_date)
df <- cbind(ids = rownames(df), df)
rownames(df) <- 1:nrow(df)
df$mecA



##create dataframe for graphing
lengths_nodes<- node.depth.edgelength(phy)[-1:-335]#calculates length of branch or year. You negative to get the nodes lengths since first node start after last tip. Number of Nodes is 1-tips so working from the back will get nodes lengths. If tips are added/removed change numbers.

#i think length nodes may be getting messed up with making nodes near zero or various factor. Probably take length as a more relative measure for now.
ancester<-ace(df$mecA,phy,type="discrete")


df_phy<-data.frame(length_node= lengths_nodes+1584,liklihood =ancester$lik.anc[,2])
#constant is added to get years. Pretty sure is safe assumption
df_phy2<-df_phy[order(df_phy$length_node),]#orders so earlier nodes are first

new_dates<-as.Date(as.yearmon(round(df_phy2$length_node, digits=0)))#just used for the runner function


#calculated aveareg for ten year window. Look at runner package if confused. 
runner_new<-runner::runner(
  x= df_phy2$liklihood,
  k = "10 years",
  idx = new_dates,
  f = function(x) mean(x)
)

q<- ggplot(df_phy2, aes(x= length_node,y= liklihood), xlab= "year", ylab= 'likelihood')+ geom_point()+geom_line(aes(y=runner_new, col="moving"))

ggplotly(q)



##how to add tip. Recomending looking at Liam Revell stuff if confused or error.

tip<-list(edge=matrix(c(2L,1L),1,2),tip.label="added_sample", edge.length=0,Nnode=1L)#edge length here is zero since trying to add trait


class(tip)<-"phylo"
phy$Nnode
phy<-bind.tree("phylo" ,tip,where=337)#add at desired node

plotTree(phy, node.numbers=T)
phy$Nnode
Ntip(xtree)

xtree<-bind.tree(xtree, tip, where=363)#bind again, but be careful node numbers change

plotTree(xtree, node.numbers=T)
xtree$tip.label

#tree@ is a useful function 

