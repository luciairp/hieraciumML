library(vegan)
library(janitor)
library(tidyverse)

datos <- read.csv("hieracium1819.csv", header=T, sep = ";")

row.names(datos) <- datos[,1]
datos <- datos[,-1]

datos <- clean_names(datos)
colnames(datos) <- abbreviate(colnames(datos),minlength=6)

# escalaBB <- c("0", "r", "+", "1", "2", "3", "4", "5")
# recodifico según 

datos[datos == "r" ] <- 0.01
datos[datos == "+" ] <- 0.1
datos[datos == 1 ] <- 2.5
datos[datos == 1 ] <- 2.5
datos[datos == 2 ] <- 17.5
datos[datos == 3 ] <- 37.5
datos[datos == 4 ] <- 62.5
datos[datos == 5 ] <- 87.5
datos[] <- lapply(datos, as.numeric)

datos.t <- t(datos)

# 1 - riqueza -------------------------------------------------------------

riqueza <- apply(datos>0,MARGIN = 1,sum) # "1" es por fila
riqueza

barplot(riqueza)
barplot(riqueza,cex.names = .7, col = 'darkgreen',las = 2)

datos.z <- read.csv("zonas_tiempos.csv", header=T, sep = ";")
datos.z$zona <- as.factor(datos.z$zona)

specpool(x = datos,pool = datos.z$zona)

# por tiempo
datos.z$tiempo <- as.factor(datos.z$tiempo)
specpool(x = datos,pool = datos.z$tiempo)

datos.z$zt <- interaction(datos.z$zona,datos.z$tiempo)

specpool(x = datos, pool = datos.z$zt)

# 2 - diversidad ----------------------------------------------------------

diversity(datos,index="shannon", MARGIN = 1) #spp x sitio
diversity(datos.t,index="shannon", MARGIN = 2) # sito x spp

datos.agg <- aggregate(datos, by = list(datos.z$zona), FUN='sum')
row.names(datos.agg) <- datos.agg[,1]
datos.agg <- datos.agg[,-1]

diversity(datos.agg,MARGIN = 1, index = "shannon")
diversity(datos.agg,MARGIN = 1, index = "simpson")

# por zona y tiempo

datos.agg.zt <- aggregate(datos, by = list(datos.z$zt), FUN='sum')
row.names(datos.agg.zt) <- datos.agg.zt[,1]
datos.agg.zt <- datos.agg.zt[,-1]

diversity(datos.agg.zt,MARGIN = 1, index = "shannon")
diversity(datos.agg.zt,MARGIN = 1, index = "simpson")


# 3 - curvas RA -----------------------------------------------------------

curva.ra <- radfit(datos.t)
plot(curva.ra)

curva.amb <- radfit(datos.agg)
plot(curva.amb)

summary(curva.amb)

curva.amb.zt <- radfit(datos.agg.zt)
plot(curva.amb.zt)

summary(curva.amb.zt)


# 4 - PCA -----------------------------------------------------------------

pca <- rda(datos)

pca
eigenvals(pca)

decostand(datos.agg,method="max",digits=2)

pca <- rda(datos, scale = T)

screeplot(pca)

# scaling 1 a sitios, mejor en este caso
plot(pca,display = "sites",type="text",scaling=1)
text(pca,display="species",labels=colnames(datos),col="red",cex=0.6)
scores(pca,choices=1:2,"sites", scaling=1)
scores(pca,choices=1:2,"species",scaling=1)

plot(pca,display = "sites",type="text",scaling=2)
text(pca,display="species",labels=colnames(datos),col="red",cex=0.6)
scores(pca,choices=1:2,"sites", scaling=2)
scores(pca,choices=1:2,"species",scaling=2)



# 5 - CCA -----------------------------------------------------------------

canon <- cca(datos,datos.z[,c(2,3)])

canon

op <- plot(canon,type="n")
text(canon,display="bp",col="blue")
text(canon,display="species",col="red",cex=0.7)
text(canon,display="sites",cex=0.8)

scores(canon,display="sites")
scores(canon,display="species")
scores(canon,display="bp")

inertcomp(canon,display="sites") 
goodness(canon,choices=1:2,statistic="explained",display="sites")
anova(canon,test="F")

# por zona y cca con tiempo
canon <- cca(datos[1:40,],datos.z[1:40,c(2,3)])
canon

op <- plot(canon,type="n")
text(canon,display="bp",col="blue")
text(canon,display="species",col="red",cex=0.7)
text(canon,display="sites",cex=0.8)

canon <- cca(datos[41:80,],datos.z[41:80,c(2,3)])
canon

op <- plot(canon,type="n")
text(canon,display="bp",col="blue")
text(canon,display="species",col="red",cex=0.7)
text(canon,display="sites",cex=0.8)

canon <- cca(datos[81:120,],datos.z[81:120,c(2,3)])
canon

op <- plot(canon,type="n")
text(canon,display="bp",col="blue")
text(canon,display="species",col="red",cex=0.7)
text(canon,display="sites",cex=0.8)


# 6 - MDS Bray-Curtis ---------------------------------------------------------------

disbc <- vegdist(datos, method = "bray", diag = T, upper = T)

mds.bc <- cmdscale(disbc, eig = TRUE, k = 2)
mds.bc$GOF

x <- mds.bc$points[,1]
y <- mds.bc$points[,2]
plot(x, y, xlab="Coordenada 1", ylab="Coordenada 2", 
     main="MDS Bray-Curtis",	type="n")
abline(h=0,v=0)

color <- c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),rep(6,20))

text(x, y, labels = datos.z$sitio, cex=.7, col = )


# 7 - NM MDS --------------------------------------------------------------

nmds <- metaMDS(datos,distance = "bray", k = 5,trymax=100)
stressplot(nmds)

ordiplot(nmds,type="n")
ordihull(nmds,groups=color,draw="polygon",col=c(1,2,3,4,5,6),label=F)
orditorp(nmds,display="species",col="blue",air=0.1,cex=1)
orditorp(nmds,display="sites",col=color,air=0.1,cex=.75)
ordiellipse(nmds,groups=color, kind = 'sd') # opciones 'ehull' convexhull, se, sd

ordiplot(nmds,type="points")

ordiplot(nmds,type="n")
nmds <- metaMDS(datos[1:40,],distance = "bray", k = 5,trymax=100)
ordihull(nmds,groups=color[1:40],draw="polygon",col=c(1,2,3,4,5,6),label=F)

ordiplot(nmds,type="n")
nmds <- metaMDS(datos[41:80,],distance = "bray", k = 5,trymax=100)
ordihull(nmds,groups=color[41:80],draw="polygon",col=c(1,2,3,4,5,6),label=F)

ordiplot(nmds,type="n")
nmds <- metaMDS(datos[81:120,],distance = "bray", k = 5,trymax=100)
ordihull(nmds,groups=color[81:120],draw="polygon",col=c(1,2,3,4,5,6),label=F)


# 8 - clasificación -------------------------------------------------------------

# disimilitud

betadiver(help=TRUE)

betadiver(datos,method="g")
betadiver(datos,method="gl") #dif dada por cant de especies -riqueza- unicas
betadiver(datos,method="j") #ojo: similitud de Jaccard
1-betadiver(datos,method="j") #disimilitud de Jaccard

# por área y tiempo
# uso datos agregados
betadiver(datos.agg.zt,method="g")
betadiver(datos.agg.zt,method="gl")
betadiver(datos.agg.zt,method="j") # similitud Jaccard

distbetad <- betadiver(datos.agg.zt)
plot(distbetad)

(euc <- dist(datos.agg.zt, method = "euclidean"))
(b <- vegdist(datos.agg.zt, method = "bray", diag = TRUE))
(gow <- vegdist(datos.agg.zt, method = "gower", diag = TRUE))
(moris <- vegdist(datos.agg.zt, method = "morisita", diag = TRUE))
(horn <- vegdist(datos.agg.zt, method = "horn", diag = TRUE))

euc.cl <- hclust(euc, method = "complete")
str(euc.cl)
plot(euc.cl)

b.cl <- hclust(b, method = "complete")
plot(b.cl)
rect.hclust(b.cl,h=0.3) #dibuja grupos con dist=h

cutree(b.cl,k=3) # si 3 grupos, cuales muestran quedan en cada uno?
cutree(b.cl,k=1:6) # cual k recupera la estructura original de grupos?

cluster.km <- kmeans(datos, centers=3, nstart = 25)
cluster.km$cluster

plot(datos, col = cluster.km$cluster)
points(cluster.km$centers, col = 1:3, pch = 8)



# 9 - indicadoras ---------------------------------------------------------

library(indicspecies)
indval <- multipatt(datos, datos.z$zt,control = how(nperm=999))
summary(indval, indvalcomp=TRUE)

# By setting alpha = 1 we want to display the group to which 
# each species is associated, regardless of whether the 
# association significant or not.

indval$sign # muestra asoc de todas las spp con grupos

# con datos de pres-aus
datos.pa <- ifelse(datos > 0,1,0)
phi <-multipatt(datos.pa, datos.z$zt, func = "r", control = how(nperm=999))
#si diferente número de muestras en grupos entonces usar func r.g.

summary(phi)
round(head(phi$str),3)

# sin combinación de grupos: duleg = T
indval.singr <- multipatt(datos, datos.z$zt,control = how(nperm=999),duleg = T)
summary(indval.singr, indvalcomp=TRUE)

# restringiendo grupos según singletons y pares: max.order = 2
indvalrest <- multipatt(datos, datos.z$zt, max.order = 2, control = how(nperm=999))
summary(indvalrest)
indvalrest$sign

# elegir grupos para comparar: restcomb = c(...) numerados
indval.sel <- multipatt(datos, datos.z$zt, restcomb = c(1:6,22,41), control = how(nperm=999))
summary(indval.sel)
indval.sel$sign

# con datos pa y eligiendo combinaciones
datos.pa <- ifelse(datos > 0,1,0)
phi <-multipatt(datos.pa, datos.z$zt, func = "r", control = how(nperm=999),restcomb = c(1:6,22,41))
#si diferente número de muestras en grupos entonces usar func r.g.

summary(phi)
round(head(phi$str),3)

# uso combinaciones de spp como indicadoras - funcion combinespecies() 
datoscomb <- combinespecies(datos, max.order = 3)$XC 
#construye objeto con sitios x spp y combinaciones de spp 
dim(datoscomb)

indvalspcomb <- multipatt(datoscomb, datos.z$zona, duleg = TRUE, control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)


# especies indicadoras por grupo
sc <- indicators(datos, cluster=datos.z$zt, group="humedo.18",max.order = 3, verbose=TRUE,At=0.5, Bt=0.2)
print(sc, sqrtIVt = 0.6)

sc <- indicators(datos, cluster=datos.z$t, group="18",max.order = 3, verbose=TRUE,At=0.5, Bt=0.2)

print(sc, sqrtIVt = 0.6)



