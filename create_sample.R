library(plyr)
library(staphopia)
library(zoo)

TOKEN = '56166f58c4617fd2c4357794b81cd806d8a3a055'
##insert token

all<-get_all_tags()
all
##look at tags and choose samples wanted

##put in sample tags if method wanted
list_of_tags<- c(218,222,1131,283,1134,267,539)
#this is a list of tags which have genetic diversity for example


sample_id_list<- vector(mode="list", length=0)
for(i in 1:length(list_of_tags)){
samples<-get_samples_by_tag(list_of_tags[i])
sample_id_list<-append(sample_id_list, samples$sample_id)
}

#will allow to get all sample_ids in one list


sample_id_list<-unlist(sample_id_list)
sample_id_list<-sample(sample_id_list,1500)
meta_large_sample <- get_metadata(sample_id_list)
sample_id_list_meta <- unique(subset(meta_large_sample, collection_date != ""))
#process data to only get samples with collection_date. This should create sample size smaller than 1500



all_genes <- get_genes(sample_id_list_meta$sample_id, exclude_sequence = T)
allgenesn <- subset(all_genes, name != 'none')
genetable <- table(allgenesn$sample_id, allgenesn$name)
genetable[which(genetable > 1)] <- 1
df<-as.data.frame.matrix(genetable) 
df$date<-as.numeric(substr(sample_id_list_meta$collection_date,1,4))
df$sequence_type<-get_sequence_type(sample_id_list_meta$sample_id)
df <- cbind(ids = rownames(df), df)
rownames(df) <- 1:nrow(df)
#####create dataframe with genes, dates and idsm and Sequence types

##basic analysis

hist(as.numeric(df$date)) #hisogram of dates


logistic <- glm(df$mecA ~df$date,  binomial(link = "logit"))
model<-predict(logistic,data.frame(df$date) ,type="response")
plot(df$date,df$mecA)
lines(df$date,model)
#logisitic model is a little off.



xtabs(~df$date+df$mecA) #gives per year number of samples that have mecA and dont. 0 is absence 1 is presence
xtabs(df$mecA ~ df$date) #gives per year number of sample with given gene



test<-aggregate(df[,1:length(colnames(df))],by=list(df$date), mean)
qplot(test$date, test$mecA)
#calculates percent of samples in given year with mecA

#note any gene can be tested

##end code with dataframe that can be used for further analysis. Not much new here
df