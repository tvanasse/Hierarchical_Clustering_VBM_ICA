# The code below reproduces clustering analyses and figures conducted in the manuscript "BrainMapVBM: 
# An Environment for Structural Meta-analysis" by Vanasse et al. (Figure 5, Sup. Figure 6, Sup. Figure 7) 
# For details regarding the motivation behind analyses and the interpretation of results, see the manuscript.
#
# written by Thomas Vanasse (thomas.j.vanasse@gmail.com), June 2016 - November 2017

library(gplots)
library(cluster)
library(dendextend)
library(RColorBrewer)

# set working directory
setwd("path to directory")

# upload matrix of network-disease lodings (43 diseases x 21 networks), mat
mat <- data.matrix(read.csv("data_mat.csv", header = FALSE, sep = ","))

# add row/column names, icd_names & col_names
icd_names <- read.csv("ICD_labels.csv", header = FALSE, sep=",")
col_names <- read.csv("component_labels.csv", header = FALSE, sep=",")
rownames(mat) <- t(icd_names)
colnames(mat) <- t(col_names)

# row distance matrix, rd
rd<-as.dist(1 - cor(t(mat), method="pearson"))

# row cluster tree, rc
rc<-hclust(rd, "average")

# column distance matrix, cd
cd<-as.dist(1 - cor(mat, method="pearson"))

# column cluster tree, cc
cc<-hclust(cd, "average")

# create row and column dendrograms, den1 & dend 2 respectively
dend1<-as.dendrogram(rc)
dend2<-as.dendrogram(cc)

# create color scale of heatmap, colors
max(mat)
min(mat)
colors <- c(seq(-0.05, 0.024, length.out=40),seq(0.025, 0.18,length.out=80))

#set color pallette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# gap statistic plots
par(mfrow=c(1,2))

# component gap statistic plots
mydist <- function(x) as.dist((1-cor(t(x))))
mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "average"),k=k))
myclusGap <- clusGap(t(mat),
                     FUN = mycluster, 
                     K.max = 18, 
                     B = 100)
plot(myclusGap, main="Component Gap Statistic", ylim=c(-0.10,0.30))
abline(v=c(9), lty=2)

# diseases gap statistic plots
mydist <- function(x) as.dist((1-cor(t(x))))
mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "average"),k=k))
myclusGap <- clusGap(mat,
                     FUN = mycluster, 
                     K.max = 30, 
                     B = 100)
plot(myclusGap, main="Disease Gap Statistic")
abline(v=c(11), lty=2)

# silhouette plots 
par(mfrow=c(1,2))
require(cluster)

# component silhouette clustering
hsilo=c()
for (i in 2:16){
  hsil <- silhouette(cutree(cc, k = i), cd)
  hsilo[i]<-summary(hsil)$avg.width
  print(hsilo[i])
}
plot(hsilo, ylab = "Average Silhouette Width",
     xlab="K Clusters", main="Components")
abline(v=c(9), lty=2)
lines(hsilo)

# component silhouette lengths, k = 9
sil_cl <- silhouette(cutree(cc, k = 9), cd)
rownames(sil_cl) <- t(col_names)
plot(sil_cl, main = "Silhouette of Component Clusters")

# disease silhouette clustering
hsilo_dis=c()
for (i in 2:35){
  hsil <- silhouette(cutree(rc, k = i), rd)
  hsilo_dis[i] <- summary(hsil)$avg.width
  print(hsilo_dis[i])
}
plot(hsilo_dis, ylab = "Average Silhouette Width", 
     xlab="K Clusters", main = "Diseases")
abline(v=c(11),lty=2)
lines(hsilo_dis)

# disease silhouette lengths, k = 11
sil_cl <- silhouette(cutree(rc, k = 11), rd)
rownames(sil_cl) <- t(icd_names)
plot(sil_cl, main = "Silhouette of Disease Clusters")

# set the colors of k branches
dend1 <- color_branches(dend1, k =11, col = col_vector)
# set colors for k branches (columns)
dend2 <- color_branches(dend2, k = 9, col = col_vector)

# set colors of "branches" in dendrogram
col_labels <- get_leaves_branches_col(dend1)
column_col_labels <- get_leaves_branches_col(dend2)
col_labels <- col_labels[order(order.dendrogram(dend1))]
column_col_labels <- column_col_labels[order(order.dendrogram(dend2))]

# plot heatmap w/ dendrograms
heatmap.2(mat, 
          Rowv=dend1,
          Colv=dend2,
          breaks=colors,
          symm=F, 
          symkey=F, 
          margins=c(5,18), 
          trace="none", 
          density.info='histogram', 
          RowSideColors=col_labels,
          ColSideColors=column_col_labels,
          cexRow=0.98, #font: cex = 0.2 + 1/log10(nr)
          cexCol=1.3) # to add nice colored strips   

# plot heatmap w/o dendrograms
heatmap.2(mat, 
          dendrogram = 'none',
          Rowv = FALSE,
          Colv = FALSE,
          breaks=colors,
          symm=F, 
          symkey=F, 
          margins=c(5,18), 
          trace="none", 
          key = FALSE,
          #density.info='histogram', 
          cexRow=0.98, #font: cex = 0.2 + 1/log10(nr)
          cexCol=1.3) 
#export 1,100 x 800
#export 1,300 x 800

# wilcoxon rank-sum test (psychiatric [0] vs. neurological disease [1] weights)
mat <- data.matrix(read.csv("psychiatric_vs_neurological_rank_sum.csv", header = FALSE, sep = ","))
wilcox.test(mat[,2]~mat[,1])

#remove graphic
dev.off()



