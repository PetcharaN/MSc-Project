#This script calculates the Pearson correlation coefficient for comparing the alternative splicing events in wheat circadian clock genes to modules defined by Rees et al

setwd("~/Project/Pearson modules")

eigengenes <- read.table('eigengene_averages.txt',header = TRUE, row.names = 1)
bin_events <- read.table('bin_events.txt',header = TRUE,row.names = 1)

rownames(eigengenes)
rownames(bin_events)

cor_stat <- matrix(data = 0, ncol=ncol(bin_events), nrow=ncol(eigengenes))
cor_pval <- matrix(data = 0, ncol=ncol(bin_events), nrow=ncol(eigengenes))


rownames(cor_stat) <- colnames(eigengenes)
rownames(cor_pval) <- colnames(eigengenes)

colnames(cor_stat) <- colnames(bin_events)
colnames(cor_pval) <- colnames(bin_events)


for (i in 1:ncol(bin_events)) 
{
  for (j in 1:ncol(eigengenes))
  {
    a <- cor.test(eigengenes[,j], as.numeric(bin_events[,i]))
    a$estimate
    a$p.value
    
    cor_stat[j,i] <- a$estimate
    cor_pval[j,i] <- a$p.value
  }
}

write.table(cor_stat, "cor_stats_WGCNA.txt", row.names=T, quote=F, col.names = T, sep="\t")
write.table(cor_pval, "cor_stats_WGCNA.txt", row.names=T, quote=F, col.names = T, sep="\t")
