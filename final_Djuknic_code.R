library(ggplot2)
library(vegan)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(clValid)
library(geosphere)

table0=read.table("final_Djuknic_for_analyses_2bioclimLU.txt", header=TRUE)

### 
# abundance and frequency and Rapoport's rule
###

sp_abundance=c()
sp_frequency=c()
lista=list()
range_lista=c()
mean_lista=c()

for (i in 11:47){
  sp_abundance[i]=sum(table0[,i])
  sp_frequency[i]=0
  elev_rang=c()
  for (j in 1:nrow(table0)){
    if(table0[j,i]>0){
      sp_frequency[i]=sp_frequency[i]+1
      elev_rang=c(elev_rang, table0$nadmorska_visina_.m.[j])
    }
  }
  lista[[i]]=elev_rang
  range_lista[i]=max(elev_rang)-min(elev_rang)
  mean_lista[i]=mean(max(elev_rang),min(elev_rang))
}

#Fig.1b
ggplot(data=NULL) + geom_point(aes(x=100*sp_frequency[10:55]/255, y=sp_abundance[10:55]), show.legend = FALSE) + geom_point(aes(x=sp_frequency[sp_frequency>25&sp_abundance>1000]*100/255, y=sp_abundance[sp_frequency>25&sp_abundance>1000]), col="#E41A1C", size=3, show.legend = FALSE) + geom_point(aes(x=sp_frequency[sp_frequency<25&sp_abundance>1000]*100/255, y=sp_abundance[sp_frequency<25&sp_abundance>1000]), col="#377EB8", size=3, show.legend = FALSE) + geom_point(aes(x=sp_frequency[sp_frequency>25&sp_abundance<1000]*100/255, y=sp_abundance[sp_frequency>25&sp_abundance<1000]), col="#4DAF4A", size=3, show.legend = FALSE) + xlab("Site frequency (%)") + ylab("Total species abundance")

###
# alpha diversity
###

loc_abundance=rowSums(table0[,11:47])
richness=specnumber(table0[,11:47])
shannons=diversity(table0[,11:47], index="shannon")
simpsons=diversity(table0[,11:47], index="invsimpson")

Renyi_div= as.data.frame(renyi(table0[,11:47]))
Renyi_div=cbind(rownames(table0), Renyi_div)

zacrtanje=melt(Renyi_div, id.vars = c("rownames(table0)"), variable.name="scale", value.name="value")
colnames(zacrtanje)[1]=c("locality")

#Fig. 1c
qplot(scale, value, data=zacrtanje, col=as.factor(rep(richness,times=11)), group=locality) + geom_line(alpha = 0.5)

#Fig. 1d
qplot(table0$E_long, table0$N_lat, size=richness, xlab="Longitude", ylab="Latitude")  + geom_point(aes(colour=table0$LU_vrednosti)) + scale_colour_gradient(low="lightyellow", high = "deeppink4")

###
# NMDS
###

ord1 <- metaMDS(table0[,11:47], distance="bray", trymax = 1000)
#check the NMDS plot
plot(ord1)

# Fig. 3a
qplot(ord1$points[,1], ord1$points[,2], col=table0$nadmorska_visina_.m.) + geom_point(size=4) + theme_classic()
# Fig. 3b
qplot(table0$nadmorska_visina_.m.[table0$Sliv_kod!=4], ord1$points[table0$Sliv_kod!=4,1], col=as.factor(table0$Sliv_kod[table0$Sliv_kod!=4])) + geom_smooth(method="lm") + theme_classic()

### 
# Mantel test
### 

#abundance data frame
abund = table0[,11:47]
#environmental vector
env = jeca1$nadmorska_visina_.m.
#longitude and latitude 
geo = data.frame(jeca1$E_long, jeca1$N_lat)
#abundance data frame - Cao dissimilarity
dist.abund = vegdist(abund, method = "cao")
#ploting abundance distances among sites
dend <- as.dendrogram(hclust(dist.abund))
dend <- color_labels(dend, h = 6, col=brewer.pal(9,"Paired"))
#Fig.7
plot(dend)
#localities' cluster identity
klasteri=cutree(dend, h=6)

cor.test(mean_lista[sp_frequency>1],range_lista[sp_frequency>1])
qplot(mean_lista[sp_frequency>1],range_lista[sp_frequency>1]) + geom_smooth(method="lm", level=0.99) + geom_text_repel(aes(label=colnames(jeca1)[sp_frequency>1]))
