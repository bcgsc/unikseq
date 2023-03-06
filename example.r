library(ggplot2)

dfa<-read.table("unikseq-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("C. maximus"), " Mt genome"))

# Stacked
ggplot(dfa, aes(fill=condition, y=value, x=position)) + 
  geom_col() + labs(x=my_x_title) + ylab("Proportion of species with 25-mers") + coord_flip()


======
PLOT ON A LOG10 SCALE

library(ggplot2)
library(ggallin)
library(scales)

dfa<-read.table("unikseq-r_CEMA.fa-i_shark.fa-o_actinopterygii.fa-k25-uniqueKmers.tsv", sep="\t", header = TRUE)
my_x_title <- expression(paste("Position of 25-mers on ", italic("C. maximus"), " Mt genome"))

# Stacked
ggplot(dfa, aes(fill=condition, y=value*1000, x=position)) + 
  scale_y_continuous(trans = pseudolog10_trans,breaks=c(-1000,-100,-10,-1,0,1,10,100,1000)) +
  geom_col() + labs(x=my_x_title) + ylab("Proportion (x10) of species with 25-mers") + coord_flip()

