setwd("/Users/thomasvanasse/Google Drive/RESEARCH/METHODS/Methods_work/SHARED_DATA/HCA/")

library(gplots)
library(cluster)

#upload matrix of network-disease lodings (21 networks)
mat <- data.matrix(read.csv("data_mat.csv", header = FALSE, sep = ","))


#Add row/column namess
icd_names <- read.csv("ICD_labels.csv", header = FALSE, sep=",")
col_names <- read.csv("component_labels.csv", header = FALSE, sep=",")
rownames(mat)<-t(icd_names)
colnames(mat)<-t(col_names)

#Row distance matrix
dissimilarity <- 1 - cor(t(mat), method="pearson")
rd <- as.dist(dissimilarity)

#rd<-dist(mat)
rc<-hclust(rd, "average")
#cophenetic distance matrix, cd
cd<-dist(t(mat))

#Column distance matrix
dissimilarity_c <- 1 - cor(mat, method="pearson")
cd <- as.dist(dissimilarity_c)

#cd<-dist(mat)
cc<-hclust(cd, "average")
cd_c<-dist(mat)

#Cophenetic Correlation 
cor(rd, cophenetic(rc))
cor(cd, cophenetic(cc))

#Create row and column dendrograms
dend1 <- as.dendrogram(rc)
dend2 <- as.dendrogram(cc)

#Color scale of heatmap
max(mat)
min(mat)
colors<- c(seq(-0.05, 0.024, length.out=40),seq(0.025, 0.18,length.out=80))

# Get the dendextend package
if(!require(dendextend)) install.packages("dendextend")
library(dendextend)
# get some colors
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Gap Statistic
library(cluster)
par(mfrow=c(1,2))
#Component
mydist <- function(x) as.dist((1-cor(t(x))))
mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "average"),k=k))
myclusGap <- clusGap(t(mat),
                     FUN = mycluster, 
                     K.max = 18, 
                     B = 100)
plot(myclusGap, main="Component Gap Statistic", ylim=c(-0.10,0.30))
abline(v=c(9), lty=2)
#Disease
mydist <- function(x) as.dist((1-cor(t(x))))
mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "average"),k=k))
myclusGap <- clusGap(mat,
                     FUN = mycluster, 
                     K.max = 30, 
                     B = 100)
plot(myclusGap, main="Disease Gap Statistic")
abline(v=c(11), lty=2)


#Silhouette
par(mfrow=c(1,2))
require(cluster)
#component silhouette clustering
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

#Silhouette Plot
sil_cl <- silhouette(cutree(cc, k = 9), cd)
rownames(sil_cl) <- t(col_names)
plot(sil_cl, main = "Silhouette of Component Clusters")

#Disease Clustering
hsilo_dis=c()
for (i in 2:35){
  hsil <- silhouette(cutree(rc, k = i), rd)
  hsilo_dis[i]<-summary(hsil)$avg.width
  print(hsilo_dis[i])
}
plot(hsilo_dis, ylab = "Average Silhouette Width", 
     xlab="K Clusters", main = "Diseases")
abline(v=c(11),lty=2)
lines(hsilo_dis)

#Silhouette Plot
c23 <- col_vector[1:11]
sil_cl <- silhouette(cutree(rc, k = 11), rd)
rownames(sil_cl) <- t(icd_names)
plot(sil_cl, main = "Silhouette of Disease Clusters")


# Set the colors of k branches
dend1 <- color_branches(dend1, k =11, col = col_vector)
# Set colors for k branches (columns)
dend2 <- color_branches(dend2, k = 9, col = col_vector)



source("https://raw.githubusercontent.com/talgalili/dendextend/master/R/attr_access.R")
col_labels <- get_leaves_branches_col(dend1)
column_col_labels <- get_leaves_branches_col(dend2)
# But due to the way heatmap.2 works - we need to fix it to be in the 
# order of the data!    
col_labels <- col_labels[order(order.dendrogram(dend1))]
column_col_labels <- column_col_labels[order(order.dendrogram(dend2))]


# Creating Heat Map
if(!require(gplots)) install.packages("gplots")
library(gplots)

heatmap.2(mat, 
          Rowv=dend1,
          Colv=dend2,
          breaks=colors,
          symm=F, 
          symkey=F, 
          margins=c(5,18), 
          trace="none", 
          density.info='histogram', 
          #na.rm=TRUE,
          RowSideColors=col_labels,
          ColSideColors=column_col_labels,
          cexRow=0.98, #font: cex = 0.2 + 1/log10(nr)
          cexCol=1.3) # to add nice colored strips        
#colRow = col_labels # to add nice colored labels - only for qplots 2.17.0 and higher

#remove graphic
dev.off()



